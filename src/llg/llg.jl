"""
GPU kernel to compute the standard LLG equation,
    dm/dt = - gamma_L * (m x H) - alpha*gamma_L* m x (m x H)
where gamma_L = gamma/(1+alpha^2).
"""
@kernel function llg_rhs_kernel!(dm_dt, @Const(m), @Const(h), @Const(pins), alpha::T,
                                 gamma::T, precession::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds pin::Bool = pins[I]

    if pin
        @inbounds dm_dt[j] = 0
        @inbounds dm_dt[j + 1] = 0
        @inbounds dm_dt[j + 2] = 0
    else
        a::T = -gamma / (1 + alpha * alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm::T = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz

        f1, f2, f3 = T(0), T(0), T(0)
        if precession
            f1 = cross_x(mx, my, mz, h1, h2, h3)
            f2 = cross_y(mx, my, mz, h1, h2, h3)
            f3 = cross_z(mx, my, mz, h1, h2, h3)
        end

        @inbounds dm_dt[j] = a * (f1 - h1 * alpha)
        @inbounds dm_dt[j + 1] = a * (f2 - h2 * alpha)
        @inbounds dm_dt[j + 2] = a * (f3 - h3 * alpha)
    end
end

"""
LLG call_back function that will be called by the integrator.
"""
function llg_call_back(sim::AbstractSim, dm_dt::AbstractArray{T,1},
                       spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    driver = sim.driver

    effective_field(sim, spin, t)

    @timeit MicroMagnetic.timer "llg" begin
        kernel! = llg_rhs_kernel!(default_backend[], groupsize[])
        kernel!(dm_dt, spin, sim.field, sim.pins, T(driver.alpha), T(driver.gamma),
                driver.precession; ndrange=N)
    end

    return nothing
end


@kernel function spatial_llg_rhs_kernel!(dm_dt, @Const(m), @Const(h), @Const(pins), @Const(alphas),
                                 gamma::T, precession::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds pin::Bool = pins[I]
    @inbounds alpha::T = alphas[I]

    if pin
        @inbounds dm_dt[j] = 0
        @inbounds dm_dt[j + 1] = 0
        @inbounds dm_dt[j + 2] = 0
    else
        @inbounds a::T = -gamma / (1 + alpha * alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm::T = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz

        f1, f2, f3 = T(0), T(0), T(0)
        if precession
            f1 = cross_x(mx, my, mz, h1, h2, h3)
            f2 = cross_y(mx, my, mz, h1, h2, h3)
            f3 = cross_z(mx, my, mz, h1, h2, h3)
        end

        @inbounds dm_dt[j] = a * (f1 - h1 * alpha)
        @inbounds dm_dt[j + 1] = a * (f2 - h2 * alpha)
        @inbounds dm_dt[j + 2] = a * (f3 - h3 * alpha)
    end
end

"""
LLG call_back function that will be called by the integrator.
"""
function spatial_llg_call_back(sim::AbstractSim, dm_dt::AbstractArray{T,1},
                       spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    driver = sim.driver

    effective_field(sim, spin, t)

    kernel! = spatial_llg_rhs_kernel!(default_backend[], groupsize[])
    kernel!(dm_dt, spin, sim.field, sim.pins, driver.alpha, T(driver.gamma),
                driver.precession; ndrange=N)

    return nothing
end


"""
compute tau = (u.nabla) m.
"""
@kernel function field_stt_kernel!(h_stt, @Const(m), @Const(ux), @Const(uy), @Const(uz),
                                   @Const(ngbs), ut::T, dx::T, dy::T,
                                   dz::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    fx::T, fy::T, fz::T = T(0), T(0), T(0)

    #x-direction
    i1::Int32 = ngbs[1, I] #we assume that i1<0 for the area with Ms=0
    i2::Int32 = ngbs[2, I]
    # i1 * i2 may overflow
    factor::T = (i1 > 0 && i2 > 0) ? 1 / (2 * dx) : 1 / dx
    i1 < 0 && (i1 = I)
    i2 < 0 && (i2 = I)
    j1 = 3 * i1 - 2
    j2 = 3 * i2 - 2
    @inbounds u = ux[I] * factor
    @inbounds fx += u * (m[j2] - m[j1])
    @inbounds fy += u * (m[j2 + 1] - m[j1 + 1])
    @inbounds fz += u * (m[j2 + 2] - m[j1 + 2])

    #y-direction
    i1 = ngbs[3, I]
    i2 = ngbs[4, I]
    factor = (i1 > 0 && i2 > 0) ? 1 / (2 * dy) : 1 / dy
    i1 < 0 && (i1 = I)
    i2 < 0 && (i2 = I)
    j1 = 3 * i1 - 2
    j2 = 3 * i2 - 2
    @inbounds u = uy[I] * factor
    @inbounds fx += u * (m[j2] - m[j1])
    @inbounds fy += u * (m[j2 + 1] - m[j1 + 1])
    @inbounds fz += u * (m[j2 + 2] - m[j1 + 2])

    #z-direction
    i1 = ngbs[5, I]
    i2 = ngbs[6, I]
    factor = (i1 > 0 && i2 > 0) ? 1 / (2 * dz) : 1 / dz
    i1 < 0 && (i1 = I)
    i2 < 0 && (i2 = I)
    j1 = 3 * i1 - 2
    j2 = 3 * i2 - 2
    @inbounds u = uz[I] * factor
    @inbounds fx += u * (m[j2] - m[j1])
    @inbounds fy += u * (m[j2 + 1] - m[j1 + 1])
    @inbounds fz += u * (m[j2 + 2] - m[j1 + 2])

    @inbounds h_stt[j] = ut * fx
    @inbounds h_stt[j + 1] = ut * fy
    @inbounds h_stt[j + 2] = ut * fz
end

"""
Wrapper for function field_stt_kernel! [to compute tau = (u.nabla) m]
"""
function compute_field_stt(h_stt, m, ux, uy, uz, ngbs, ut::T, dx::T, dy::T, dz::T,
                           N::Int64) where {T<:AbstractFloat}
    kernel! = field_stt_kernel!(default_backend[], groupsize[])
    kernel!(h_stt, m, ux, uy, uz, ngbs, ut, dx, dy, dz; ndrange=N)

    return nothing
end

"""
GPU kernel to compute the STT torque and add it to dm/dt:
    dm/dt += (1+alpha*beta)/(1+alpha^2)*tau - (beta-alpha)/(1+alpha^2)*(m x tau)
where tau = (u.nabla) m, here tau is represented by h_stt.
"""
@kernel function llg_stt_rhs_kernel!(dm_dt, @Const(m), @Const(h), @Const(h_stt),
                                     @Const(pins), gamma::T, alpha::T,
                                     beta::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds if pins[I]
        @inbounds dm_dt[j] = 0
        @inbounds dm_dt[j + 1] = 0
        @inbounds dm_dt[j + 2] = 0
    else
        a::T = -gamma / (1 + alpha * alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm::T = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz

        f1 = cross_x(mx, my, mz, h1, h2, h3)
        f2 = cross_y(mx, my, mz, h1, h2, h3)
        f3 = cross_z(mx, my, mz, h1, h2, h3)

        @inbounds dm_dt[j] = a * (f1 - h1 * alpha)
        @inbounds dm_dt[j + 1] = a * (f2 - h2 * alpha)
        @inbounds dm_dt[j + 2] = a * (f3 - h3 * alpha)

        #the above part is the standard LLG equation

        c1::T = (1 + alpha * beta) / (1 + alpha * alpha)
        c2::T = (beta - alpha) / (1 + alpha * alpha)

        @inbounds mht = mx * h_stt[j] + my * h_stt[j + 1] + mz * h_stt[j + 2]
        @inbounds hp1 = mm * h_stt[j] - mht * mx
        @inbounds hp2 = mm * h_stt[j + 1] - mht * my
        @inbounds hp3 = mm * h_stt[j + 2] - mht * mz

        mth1 = cross_x(mx, my, mz, hp1, hp2, hp3)
        mth2 = cross_y(mx, my, mz, hp1, hp2, hp3)
        mth3 = cross_z(mx, my, mz, hp1, hp2, hp3)

        @inbounds dm_dt[j] += c1 * hp1 - c2 * mth1
        @inbounds dm_dt[j + 1] += c1 * hp2 - c2 * mth2
        @inbounds dm_dt[j + 2] += c1 * hp3 - c2 * mth3
    end
end

"""
LLG with STT extension call_back function that will be called by the integrator.
"""
function llg_stt_call_back(sim::AbstractSim, dm_dt, spin, t::Float64)
    driver = sim.driver
    mesh = sim.mesh

    effective_field(sim, spin, t)

    T = Float[]

    ut = T(driver.ufun(t))
    compute_field_stt(driver.h_stt, spin, driver.ux, driver.uy, driver.uz, mesh.ngbs, ut,
                      T(mesh.dx), T(mesh.dy), T(mesh.dz), sim.n_total)

    kernel! = llg_stt_rhs_kernel!(default_backend[], groupsize[])
    kernel!(dm_dt, spin, sim.field, driver.h_stt, sim.pins, T(driver.gamma),
            T(driver.alpha), T(driver.beta); ndrange=sim.n_total)

    return nothing
end

"""
GPU kernel to compute the STT torque for CPP case and add it to dm/dt:
    dm/dt += (a_J+alpha*b_J)/(1+alpha^2)*p_perp - (b_J-alpha*a_J)/(1+alpha^2)*(m x p_perp)
where p_perp = p - (m.p)m
"""
@kernel function llg_cpp_rhs_kernel!(dm_dt, @Const(m), @Const(h), @Const(pins), @Const(a_J),
                                     bj::T, ut::T, gamma::T, alpha::T, px::T, py::T,
                                     pz::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds if pins[I]
        @inbounds dm_dt[j] = 0
        @inbounds dm_dt[j + 1] = 0
        @inbounds dm_dt[j + 2] = 0
    else
        a::T = -gamma / (1 + alpha * alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm::T = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz

        f1 = cross_x(mx, my, mz, h1, h2, h3)
        f2 = cross_y(mx, my, mz, h1, h2, h3)
        f3 = cross_z(mx, my, mz, h1, h2, h3)

        @inbounds dm_dt[j] = a * (f1 - h1 * alpha)
        @inbounds dm_dt[j + 1] = a * (f2 - h2 * alpha)
        @inbounds dm_dt[j + 2] = a * (f3 - h3 * alpha)

        #the above part is the standard LLG equation

        @inbounds aj = ut * a_J[I]
        bj = ut * bj
        c1::T = (aj + alpha * bj) / (1 + alpha * alpha)
        c2::T = (bj - alpha * aj) / (1 + alpha * alpha)

        mpt::T = mx * px + my * py + mz * pz
        @inbounds hp1 = mm * px - mpt * mx
        @inbounds hp2 = mm * py - mpt * my
        @inbounds hp3 = mm * pz - mpt * mz

        mth1 = cross_x(mx, my, mz, hp1, hp2, hp3)
        mth2 = cross_y(mx, my, mz, hp1, hp2, hp3)
        mth3 = cross_z(mx, my, mz, hp1, hp2, hp3)

        @inbounds dm_dt[j] += c1 * hp1 - c2 * mth1
        @inbounds dm_dt[j + 1] += c1 * hp2 - c2 * mth2
        @inbounds dm_dt[j + 2] += c1 * hp3 - c2 * mth3
    end
end

"""
The call_back function for STT torque for CPP case that will be called by the integrator.
"""
function llg_cpp_call_back(sim::AbstractSim, dm_dt, spin, t::Float64)
    driver = sim.driver

    effective_field(sim, spin, t)

    T = Float[]
    ut = driver.ufun(t)
    px, py, pz = T(driver.p[1]), T(driver.p[2]), T(driver.p[3])

    kernel! = llg_cpp_rhs_kernel!(default_backend[], groupsize[])
    kernel!(dm_dt, spin, sim.field, sim.pins, driver.aj, T(driver.bj), T(ut),
            T(driver.gamma), T(driver.alpha), px, py, pz; ndrange=sim.n_total)

    return nothing
end
