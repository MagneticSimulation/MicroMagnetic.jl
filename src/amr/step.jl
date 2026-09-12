# Composite-grid time integration.
#
# Every region (base level + patches) advances with the SAME time step using
# the Gauss-Seidel projection scheme (Scheme B of García-Cervera & E,
# IEEE Trans. Magn. 39, 1766 (2003); the package's GPSM integrator):
#
#   for k = 1,2,3 (components):
#       m_k* = m_k - (m_{k+1} g_{k+2} - m_{k+2} g_{k+1})
#              - alpha (m.g) m_k + alpha |m|^2 g_k
#       hs   = anisotropy(m*) + zeeman + Phi          (pointwise + frozen demag)
#       g_k  = (I - dtg * L)^-1 (m_k* + dtg * hs_k)   (sparse Cholesky)
#
# with dtg = dt * gamma / (1 + alpha^2) and L the exchange Laplacian of the
# region (built from its ghosted stencil, so patches see the interpolated
# ghost cells as Dirichlet-like data). The magnetization is renormalized at
# the end of the step. Because the implicit operator decouples the exchange
# from the step size, the scheme is unconditionally stable and all levels
# share one dt (Section II of the 2006 paper).
#
# The stray field Phi is frozen at the beginning of the step
# (`update_demag_boxes!`); anisotropy and Zeeman are re-evaluated at every
# intermediate state, as in the package's GPSM integrator.

"""Pointwise GPSM intermediate state for component `k` (in place on `m`)."""
@kernel function amr_gpsm_star_kernel!(m, @Const(g1), @Const(g2), @Const(g3),
                                       alpha::T, k::Int) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2
    @inbounds m1 = m[j]
    @inbounds m2 = m[j + 1]
    @inbounds m3 = m[j + 2]
    @inbounds ga = g1[I]
    @inbounds gb = g2[I]
    @inbounds gc = g3[I]
    mg = m1 * ga + m2 * gb + m3 * gc
    mm = m1 * m1 + m2 * m2 + m3 * m3
    if k == 1
        @inbounds m[j] = m1 - (m2 * gc - m3 * gb) - alpha * mg * m1 + alpha * mm * ga
    elseif k == 2
        @inbounds m[j + 1] = m2 - (m3 * ga - m1 * gc) - alpha * mg * m2 + alpha * mm * gb
    else
        @inbounds m[j + 2] = m3 - (m1 * gb - m2 * ga) - alpha * mg * m3 + alpha * mm * gc
    end
end

"""
Assemble the GPSM right-hand side of component `k` for one cell of a region
box. Authoritative cells:  rhs = m_k + dtg*(hs_anis + hs_zeem + Phi)  plus the
Dirichlet injection  dtg * L_ij * g_j  from non-authoritative (Dirichlet)
neighbors, whose auxiliary field is prescribed by the composite box Gc.
Dirichlet cells themselves: rhs = Gc (the solve returns it as identity rows).
The Phi/Gc boxes live at the region's own level; box cell (a,b,c) maps to the
global level cell (i0-1+a, j0-1+b, k0-1+c) (identity for the base region),
clamped into the box for out-of-sample ghost cells.
"""
@kernel function amr_rhs_kernel!(rhs, @Const(m), @Const(hs_a), @Const(hs_z),
                                 @Const(phi), @Const(gc), @Const(isdir),
                                 k::Int, dtg::T,
                                 cx::T, cy::T, cz::T,
                                 i0::Int, j0::Int, k0::Int,
                                 ngx::Int, ngy::Int, ngz::Int,
                                 nlx::Int, nly::Int, nlz::Int) where {T<:AbstractFloat}
    a, b, c = @index(Global, NTuple)
    gi = clamp(i0 - 1 + a, 1, nlx)
    gj = clamp(j0 - 1 + b, 1, nly)
    gk = clamp(k0 - 1 + c, 1, nlz)
    Ib = _cell_index(a, b, c, ngx, ngy)
    Ig = _cell_index(gi, gj, gk, nlx, nly)
    j = 3 * Ib - 2 + (k - 1)
    e = 3 * Ig - 2 + (k - 1)
    @inbounds if isdir[Ib]
        rhs[Ib] = gc[e]
    else
        v = m[j] + dtg * (hs_a[j] + hs_z[j] + phi[e])
        # inject the prescribed g of the Dirichlet stencil neighbors
        @inbounds begin
            if a > 1 && isdir[Ib - 1]
                gj1 = clamp(i0 - 2 + a, 1, nlx)
                v += cx * gc[3 * _cell_index(gj1, gj, gk, nlx, nly) - 2 + (k - 1)]
            end
            if a < ngx && isdir[Ib + 1]
                gj1 = clamp(i0 + a, 1, nlx)
                v += cx * gc[3 * _cell_index(gj1, gj, gk, nlx, nly) - 2 + (k - 1)]
            end
            if b > 1 && isdir[Ib - ngx]
                gj2 = clamp(j0 - 2 + b, 1, nly)
                v += cy * gc[3 * _cell_index(gi, gj2, gk, nlx, nly) - 2 + (k - 1)]
            end
            if b < ngy && isdir[Ib + ngx]
                gj2 = clamp(j0 + b, 1, nly)
                v += cy * gc[3 * _cell_index(gi, gj2, gk, nlx, nly) - 2 + (k - 1)]
            end
            if c > 1 && isdir[Ib - ngx * ngy]
                gj3 = clamp(k0 - 2 + c, 1, nlz)
                v += cz * gc[3 * _cell_index(gi, gj, gj3, nlx, nly) - 2 + (k - 1)]
            end
            if c < ngz && isdir[Ib + ngx * ngy]
                gj3 = clamp(k0 + c, 1, nlz)
                v += cz * gc[3 * _cell_index(gi, gj, gj3, nlx, nly) - 2 + (k - 1)]
            end
        end
        rhs[Ib] = v
    end
end

"""Advance one region by one GPSM step; returns max |dm| over its
authoritative cells (for the max dm/dt stopping criterion)."""
function _region_step!(amr::AMRSim{T}, r::AMRRegion{T}, dtg::Float64) where {T<:AbstractFloat}
    sim = r.sim
    N = sim.n_total
    m = sim.spin
    back = get_backend(m)
    copyto!(r.prespin, m)

    # level geometry for the Phi lookup (from the region's geo cache: ghost
    # widths on refined dims only, box origin/dims, level dims)
    geo = r.geo
    i0, j0, k0 = geo.i0, geo.j0, geo.k0
    ngx, ngy, ngz = geo.ngx, geo.ngy, geo.ngz
    nlx, nly, nlz = geo.nlx, geo.nly, geo.nlz
    phi = amr.Phi[r.level]

    for k in 1:3
        kernel! = amr_gpsm_star_kernel!(back, groupsize[])
        kernel!(m, r.g1, r.g2, r.g3, amr.alpha, k; ndrange=N)

        # pointwise fields at the intermediate state
        if amr.Ku != 0
            T_vol = T(sim.mesh.volume)
            nax = Fill(T(amr.axis[1]), N)
            nay = Fill(T(amr.axis[2]), N)
            naz = Fill(T(amr.axis[3]), N)
            ku = Fill(amr.Ku, N)
            ak! = anisotropy_kernel!(back, groupsize[])
            ak!(m, r.hs_a, sim.energy, ku, nax, nay, naz, sim.mu0_Ms, T_vol; ndrange=N)
        end
        if amr.H0 !== nothing
            hx = Fill(T(amr.H0[1]), N)
            hy = Fill(T(amr.H0[2]), N)
            hz = Fill(T(amr.H0[3]), N)
            zk! = zeeman_field_kernel!(back, groupsize[])
            zk!(m, r.hs_z, sim.energy, sim.mu0_Ms, hx, hy, hz, T(sim.mesh.volume),
                T(1), T(1), T(1); ndrange=N)
        end

        rk! = amr_rhs_kernel!(back, groupsize[])
        rk!(r.rhs, m, r.hs_a, r.hs_z, phi, amr.Gc[r.level], r.isdir, k, T(dtg),
            r.cx, r.cy, r.cz, i0, j0, k0, ngx, ngy, ngz, nlx, nly, nlz;
            ndrange=(ngx, ngy, ngz))

        g = k == 1 ? r.g1 : k == 2 ? r.g2 : r.g3
        copyto!(g, r.G \ r.rhs)
    end
    normalise(m, N)

    # max |dm| over authoritative cells (host reduction over the bookkeeping
    # traversal; cells under finer patches are excluded)
    compute_dm!(r.dm, m, r.prespin, N)
    dm = Array(r.dm)
    maxdm = 0.0
    for_each_authoritative_cell(amr, r) do _, Ib, _
        dm[Ib] > maxdm && (maxdm = dm[Ib])
    end
    return maxdm
end

"""
Initialize the GPSM auxiliary fields g of every region from the current
magnetization:  g_k = (I - dtg L)^-1 (m_k + dtg hs_k),  matching the
initialization of the package's GPSM integrator. Called once after
`init_m0` and again after every remesh (new regions start from a consistent
state).
"""
function _init_gpsm!(amr::AMRSim{T}, dt::Real) where {T<:AbstractFloat}
    dtg = dt * amr.gamma / (1 + amr.alpha^2)
    for r in (amr.base, (p for ps in amr.patches for p in ps)...)
        sim = r.sim
        N = sim.n_total
        m = sim.spin
        back = get_backend(m)
        geo = r.geo
        i0, j0, k0 = geo.i0, geo.j0, geo.k0
        ngx, ngy, ngz = geo.ngx, geo.ngy, geo.ngz
        nlx, nly, nlz = geo.nlx, geo.nly, geo.nlz
        phi = amr.Phi[r.level]
        for k in 1:3
            if amr.Ku != 0
                ku = Fill(amr.Ku, N)
                ak! = anisotropy_kernel!(back, groupsize[])
                ak!(m, r.hs_a, sim.energy, ku,
                    Fill(T(amr.axis[1]), N), Fill(T(amr.axis[2]), N),
                    Fill(T(amr.axis[3]), N), sim.mu0_Ms,
                    T(sim.mesh.volume); ndrange=N)
            end
            if amr.H0 !== nothing
                zk! = zeeman_field_kernel!(back, groupsize[])
                zk!(m, r.hs_z, sim.energy, sim.mu0_Ms,
                    Fill(T(amr.H0[1]), N), Fill(T(amr.H0[2]), N),
                    Fill(T(amr.H0[3]), N), T(sim.mesh.volume),
                    T(1), T(1), T(1); ndrange=N)
            end
            rk! = amr_rhs_kernel!(back, groupsize[])
            rk!(r.rhs, m, r.hs_a, r.hs_z, phi, amr.Gc[r.level], r.isdir, k, T(dtg),
                r.cx, r.cy, r.cz, i0, j0, k0, ngx, ngy, ngz, nlx, nly, nlz;
                ndrange=(ngx, ngy, ngz))
            g = k == 1 ? r.g1 : k == 2 ? r.g2 : r.g3
            copyto!(g, r.G \ r.rhs)
        end
    end
    amr.g_initialized = true
    return nothing
end

"""
    amr_step!(amr, dt = amr.dt)

Advance the composite grid by one time step: synchronize the composite
boxes, remesh if due, rebuild the layered demag boxes (frozen over the step),
refresh the ghost cells, and advance every region with the same GPSM step.
"""
function amr_step!(amr::AMRSim, dt::Real = amr.dt)
    amr.ecache = nothing   # the saved energy row is accounted after this step
    sync_composite!(amr)
    if amr.nsteps > 0 && amr.nsteps % amr.remesh_interval == 0
        remesh!(amr)
    end
    update_demag_boxes!(amr)
    fill_ghosts!(amr)
    sync_g_boxes!(amr)
    if !amr.g_initialized
        _init_gpsm!(amr, dt)
    end

    dtg = dt * amr.gamma / (1 + amr.alpha^2)
    maxdm = 0.0
    for r in (amr.base, (p for ps in amr.patches for p in ps)...)
        maxdm = max(maxdm, _region_step!(amr, r, dtg))
    end

    amr.nsteps += 1
    amr.time += dt
    amr.maxdmdt = maxdm / dt
    amr.saver.t = amr.time
    amr.saver.nsteps = amr.nsteps
    if amr.save_data
        write_data(amr)
    end
end

"""
    relax(amr::AMRSim; maxsteps=10000, stopping_dmdt=0.1, dt=amr.dt,
          save_m_every=-1, save_m_path="m", verbose=true)

Relax the composite system with fixed time steps `dt` until the maximum
|dm/dt| over the composite grid drops below `stopping_dmdt` or `maxsteps`
steps are taken. This mirrors `relax(::AbstractSim; ...)` of the uniform
driver. `save_m_every > 0` saves the magnetization (interpolated to the
finest uniform grid) every that many steps.
"""
function relax(amr::AMRSim; maxsteps::Int=10000, stopping_dmdt::Real=0.1,
               dt::Real=amr.dt, save_m_every::Int=-1, save_m_path::String="m",
               verbose::Bool=true)
    for i in 1:maxsteps
        amr_step!(amr, dt)
        if verbose && (amr.nsteps % 100 == 0 || amr.maxdmdt < stopping_dmdt)
            @info "step=$(amr.nsteps) t=$(amr.time) max_dmdt=$(amr.maxdmdt) " *
                  "patches=$(sum(length.(amr.patches)))"
        end
        if save_m_every > 0 && amr.nsteps % save_m_every == 0
            save_ovf(amr, @sprintf("%s_%06d", save_m_path, amr.nsteps))
        end
        amr.maxdmdt < stopping_dmdt && break
    end
    return amr.maxdmdt
end

"""
    run_sim(amr::AMRSim; steps=100, dt=amr.dt, save_m_every=-1, save_m_path="m")

Run `steps` fixed-dt time steps of the composite system.
"""
function run_sim(amr::AMRSim; steps::Int=100, dt::Real=amr.dt, save_m_every::Int=-1,
                 save_m_path::String="m")
    for i in 1:steps
        amr_step!(amr, dt)
        if save_m_every > 0 && amr.nsteps % save_m_every == 0
            save_ovf(amr, @sprintf("%s_%06d", save_m_path, amr.nsteps))
        end
    end
end
