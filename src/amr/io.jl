# Energy accounting, data saver and output for AMRSim.

"""
    average_m(amr::AMRSim) -> Tuple

Volume-averaged magnetization of the composite system,
`<m> = ∫ m dV / ∫ dV`, summed over the authoritative cells of every level
(each level contributes with its own cell volume; cells under finer patches
are excluded).
"""
function average_m(amr::AMRSim{T}) where {T<:AbstractFloat}
    smx = T(0)
    smy = T(0)
    smz = T(0)
    sv = T(0)
    for r in (amr.base, (p for ps in amr.patches for p in ps)...)
        m = Array(r.sim.spin)
        nxi, nyi, nzi = length(r.ir), length(r.jr), length(r.kr)
        gx, gy, gz = r.level == 1 ? (0, 0, 0) : _ghosts(amr)
        nbx, nby = r.sim.mesh.nx, r.sim.mesh.ny
        lbx, lby = amr.dims[r.level][1], amr.dims[r.level][2]
        V = T(_level_volume(amr, r.level))
        if r.level == 1
            for c in 1:nzi, b in 1:nyi, a in 1:nxi
                Ig = _cell_index(a, b, c, lbx, lby)
                amr.covered[1][Ig] && continue
                j = 3 * _cell_index(a, b, c, nbx, nby) - 2
                smx += m[j] * V
                smy += m[j + 1] * V
                smz += m[j + 2] * V
                sv += V
            end
        else
            cov = amr.covered[r.level]
            for c in 1:nzi, b in 1:nyi, a in 1:nxi
                Ig = _cell_index(first(r.ir) - 1 + a, first(r.jr) - 1 + b,
                                 first(r.kr) - 1 + c, lbx, lby)
                cov[Ig] && continue
                j = 3 * _cell_index(a + gx, b + gy, c + gz, nbx, nby) - 2
                smx += m[j] * V
                smy += m[j + 1] * V
                smz += m[j + 2] * V
                sv += V
            end
        end
    end
    return (smx / sv, smy / sv, smz / sv)
end

"""
    init_m0(amr::AMRSim, m0; norm=true)

Initialize the magnetization of the base level from a tuple, array or
function `(i,j,k,dx,dy,dz) -> (mx,my,mz)` (same conventions as
`init_m0(::AbstractSim, ...)`), then generate the initial composite grid and
GPSM contexts. For a function `m0`, the newly created fine patches are
re-evaluated from it at their own resolution (indices are global level-ℓ
cell indices, matching the base-grid convention), so the initial state
carries the intended fine structure instead of interpolated values.
"""
function init_m0(amr::AMRSim, m0::TupleOrArrayOrFunction; norm::Bool=true)
    init_m0(amr.base.sim, m0; norm=norm)
    amr.time = 0.0
    amr.nsteps = 0
    amr.maxdmdt = 0.0
    amr.ecache = nothing
    amr.saver.header_saved = false
    amr.saver.t = 0.0
    amr.saver.nsteps = 0
    sync_composite!(amr)
    remesh!(amr)
    if isa(m0, Function)
        _reinit_patches_m0!(amr, m0, norm)
        sync_composite!(amr)
        remesh!(amr)   # re-flag on the real fine structure: init_m0 returns the final initial composite grid
    end
    return true
end

"""Evaluate an analytic `m0` on every patch interior at its own resolution
(fine structure instead of the interpolated first-guess values)."""
function _reinit_patches_m0!(amr::AMRSim, m0::Function, norm::Bool)
    for l in 2:amr.levels, p in amr.patches[l - 1]
        nxp, nyp, nzp = p.sim.mesh.nx, p.sim.mesh.ny, p.sim.mesh.nz
        gx, gy, gz = _ghosts(amr)
        hx, hy, hz = _level_h(amr.base_mesh, amr.refined, l)
        spin = Array(p.sim.spin)
        for kk in (gz+1):(nzp-gz), jj in (gy+1):(nyp-gy), ii in (gx+1):(nxp-gx)
            gi = first(p.ir) - 1 + (ii - gx)
            gj = first(p.jr) - 1 + (jj - gy)
            gk = first(p.kr) - 1 + (kk - gz)
            m = m0(gi, gj, gk, hx, hy, hz)
            I = _cell_index(ii, jj, kk, nxp, nyp)
            spin[3I-2] = m[1]
            spin[3I-1] = m[2]
            spin[3I] = m[3]
            if norm
                nm = sqrt(m[1]^2 + m[2]^2 + m[3]^2)
                if nm > 0
                    spin[3I-2] /= nm
                    spin[3I-1] /= nm
                    spin[3I] /= nm
                end
            end
        end
        copyto!(p.sim.spin, spin)
    end
    return nothing
end

"""Per-level cell volume (refined dims halve)."""
_level_volume(amr::AMRSim, l::Int) =
    prod(_level_h(amr.base_mesh, amr.refined, l))

# masked sum helpers over a region: iterate the interior, skip cells under
# finer patches
function _region_energy_sum(amr::AMRSim, r::AMRRegion, buf::AbstractArray{<:Any,1})
    nxi, nyi, nzi = length(r.ir), length(r.jr), length(r.kr)
    gx, gy, gz = r.level == 1 ? (0, 0, 0) : _ghosts(amr)
    nbx, nby = r.sim.mesh.nx, r.sim.mesh.ny
    lbx, lby = amr.dims[r.level][1], amr.dims[r.level][2]
    host = Array(buf)
    s = 0.0
    if r.level == 1
        for c in 1:nzi, b in 1:nyi, a in 1:nxi
            Ig = _cell_index(a, b, c, lbx, lby)
            amr.covered[1][Ig] && continue
            s += host[_cell_index(a, b, c, nbx, nby)]
        end
    else
        cov = amr.covered[r.level]
        for c in 1:nzi, b in 1:nyi, a in 1:nxi
            Ig = _cell_index(first(r.ir) - 1 + a, first(r.jr) - 1 + b,
                             first(r.kr) - 1 + c, lbx, lby)
            cov[Ig] && continue
            s += host[_cell_index(a + gx, b + gy, c + gz, nbx, nby)]
        end
    end
    return s
end

function _region_demag_energy(amr::AMRSim{T}, r::AMRRegion{T}) where {T<:AbstractFloat}
    phi = Array(amr.Phi[r.level])
    m = Array(r.sim.spin)
    nxi, nyi, nzi = length(r.ir), length(r.jr), length(r.kr)
    gx, gy, gz = r.level == 1 ? (0, 0, 0) : _ghosts(amr)
    nbx, nby = r.sim.mesh.nx, r.sim.mesh.ny
    lbx, lby = amr.dims[r.level][1], amr.dims[r.level][2]
    mu0Ms = Array(r.sim.mu0_Ms)
    V = T(_level_volume(amr, r.level))
    s = T(0)
    if r.level == 1
        for c in 1:nzi, b in 1:nyi, a in 1:nxi
            Ig = _cell_index(a, b, c, lbx, lby)
            amr.covered[1][Ig] && continue
            I = _cell_index(a, b, c, nbx, nby)
            j = 3 * I - 2
            e = 3 * Ig - 2
            s -= 0.5 * mu0Ms[I] * V *
                 (m[j] * phi[e] + m[j + 1] * phi[e + 1] + m[j + 2] * phi[e + 2])
        end
    else
        cov = amr.covered[r.level]
        for c in 1:nzi, b in 1:nyi, a in 1:nxi
            Ig = _cell_index(first(r.ir) - 1 + a, first(r.jr) - 1 + b,
                             first(r.kr) - 1 + c, lbx, lby)
            cov[Ig] && continue
            I = _cell_index(a + gx, b + gy, c + gz, nbx, nby)
            j = 3 * I - 2
            e = 3 * Ig - 2
            s -= 0.5 * mu0Ms[I] * V *
                 (m[j] * phi[e] + m[j + 1] * phi[e + 1] + m[j + 2] * phi[e + 2])
        end
    end
    return s
end

"""
    amr_energies(amr) -> NamedTuple

Exchange, demag, anisotropy and Zeeman energies of the composite system
(J), summed over the authoritative cells of every level. The exchange
energy uses the package's ngbs stencil per region; the demag energy is
-0.5 mu0 Ms m.Phi as in `collect_h_energy`.
"""
function amr_energies(amr::AMRSim{T}) where {T<:AbstractFloat}
    E_exch = 0.0
    E_demag = 0.0
    for r in (amr.base, (p for ps in amr.patches for p in ps)...)
        exch = r.sim.interactions[findfirst(x -> isa(x, Exchange),
                                            r.sim.interactions)]
        effective_field(exch, r.sim, r.sim.spin, 0.0)
        E_exch += _region_energy_sum(amr, r, exch.energy)
        amr.demag && (E_demag += _region_demag_energy(amr, r))
    end
    E_anis = 0.0
    E_zeem = 0.0
    if amr.Ku != 0 || amr.H0 !== nothing
        for r in (amr.base, (p for ps in amr.patches for p in ps)...)
            N = r.sim.n_total
            back = get_backend(r.sim.spin)
            if amr.Ku != 0
                ku = Fill(amr.Ku, N)
                ak! = anisotropy_kernel!(back, groupsize[])
                ak!(r.sim.spin, r.hs_a, r.sim.energy, ku,
                    Fill(T(amr.axis[1]), N), Fill(T(amr.axis[2]), N),
                    Fill(T(amr.axis[3]), N), r.sim.mu0_Ms,
                    T(r.sim.mesh.volume); ndrange=N)
                E_anis += _region_energy_sum(amr, r, r.sim.energy)
            end
            if amr.H0 !== nothing
                zk! = zeeman_field_kernel!(back, groupsize[])
                zk!(r.sim.spin, r.hs_z, r.sim.energy, r.sim.mu0_Ms,
                    Fill(T(amr.H0[1]), N), Fill(T(amr.H0[2]), N),
                    Fill(T(amr.H0[3]), N), T(r.sim.mesh.volume),
                    T(1), T(1), T(1); ndrange=N)
                E_zeem += _region_energy_sum(amr, r, r.sim.energy)
            end
        end
    end
    return (exch=E_exch, demag=E_demag, anis=E_anis, zeeman=E_zeem,
            total=E_exch + E_demag + E_anis + E_zeem)
end

"""
Cached `amr_energies` for the data saver: computed once per saved row and
shared by all five energy columns (invalidated by `amr_step!`, `remesh!` and
`init_m0`). The public `amr_energies` stays a pure re-evaluation.
"""
function _amr_energies_cached(amr::AMRSim)
    amr.ecache === nothing && (amr.ecache = amr_energies(amr))
    return amr.ecache
end

function _init_amr_saver!(amr::AMRSim)
    saver = DataSaver(amr.name * ".txt", false, 0.0, 0, [])
    push!(saver.items, SaverItem("step", "<unitless>", o -> o.nsteps))
    push!(saver.items, SaverItem("time", "<s>", o -> o.time))
    push!(saver.items, SaverItem("E_total", "<J>", o -> _amr_energies_cached(o).total))
    push!(saver.items, SaverItem("E_exch", "<J>", o -> _amr_energies_cached(o).exch))
    push!(saver.items, SaverItem("E_demag", "<J>", o -> _amr_energies_cached(o).demag))
    push!(saver.items, SaverItem("E_anis", "<J>", o -> _amr_energies_cached(o).anis))
    push!(saver.items, SaverItem("max_dmdt", "<1/s>", o -> o.maxdmdt))
    push!(saver.items, SaverItem("n_patches", "<unitless>",
                                 o -> sum(length.(o.patches))))
    push!(saver.items, SaverItem("ncells", "<unitless>", composite_ncells))
    push!(saver.items, SaverItem("E_zeeman", "<J>", o -> _amr_energies_cached(o).zeeman))
    amr.saver = saver
    return saver
end

"""Append one row to the AMR data table (mirrors `write_data(::AbstractSim)`)."""
function write_data(amr::AMRSim)
    saver = amr.saver
    if !saver.header_saved
        io = open(saver.name, "w")
        write(io, "#")
        for item in saver.items
            write(io, formatstring(item.name))
        end
        write(io, "\n#")
        for item in saver.items
            write(io, formatstring(item.unit))
        end
        write(io, "\n")
        saver.header_saved = true
        close(io)
    end
    io = open(saver.name, "a")
    for item in saver.items
        write(io, formatstring(item.result(amr)))
    end
    write(io, "\n")
    return close(io)
end

"""
    save_ovf(amr::AMRSim, fname; type=Float64)

Save the magnetization interpolated to the finest equivalent uniform grid
(the level-`levels` composite box: patch data where authoritative, the
progressive coarse interpolation elsewhere -- the visualization procedure of
Section IV of the paper).
"""
function save_ovf(amr::AMRSim, fname::String; type::DataType=Float64)
    nx, ny, nz = amr.dims[amr.levels]
    hx, hy, hz = _level_h(amr.base_mesh, amr.refined, amr.levels)
    data = Array(amr.C[amr.levels])
    return mag2ovf(data, nx, ny, nz; dx=hx, dy=hy, dz=hz, fname=fname, type=type)
end
