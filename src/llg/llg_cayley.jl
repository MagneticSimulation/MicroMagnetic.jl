
@kernel function llg_rhs_cayley_kernel!(dw_dt, @Const(m), @Const(h), @Const(omega),
                                        @Const(pins), alpha::T, gamma::T,
                                        precession::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    if pins[I]
        @inbounds dw_dt[j] = 0
        @inbounds dw_dt[j + 1] = 0
        @inbounds dw_dt[j + 2] = 0
    else
        a::T = -gamma / (1 + alpha * alpha)
        b::T = alpha * a
        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm::T = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz
        f1 = -a * h1 * precession - b * cross_x(mx, my, mz, h1, h2, h3)
        f2 = -a * h2 * precession - b * cross_y(mx, my, mz, h1, h2, h3)
        f3 = -a * h3 * precession - b * cross_z(mx, my, mz, h1, h2, h3)

        @inbounds wx = omega[j]
        @inbounds wy = omega[j + 1]
        @inbounds wz = omega[j + 2]

        wf = wx * f1 + wy * f2 + wz * f3
        @inbounds dw_dt[j] = f1 - 0.5 * cross_x(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wx
        @inbounds dw_dt[j + 1] = f2 - 0.5 * cross_y(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wy
        @inbounds dw_dt[j + 2] = f3 - 0.5 * cross_z(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wz
    end
end

"""
llg_cayley_call_back function that will be called by the integrator.
"""
function llg_cayley_call_back(sim::AbstractSim, dw_dt::AbstractArray{T,1}, t::Float64,
                              omega::AbstractArray{T,1}) where {T<:AbstractFloat}
    driver = sim.driver
    N = sim.n_total

    omega_to_spin(omega, sim.prespin, sim.spin, N)
    effective_field(sim, sim.spin, t)

    kernel! = llg_rhs_cayley_kernel!(default_backend[], groupsize[])
    kernel!(dw_dt, sim.spin, sim.field, omega, sim.pins, driver.alpha, driver.gamma,
            driver.precession; ndrange=N)
    KernelAbstractions.synchronize(default_backend[])

    return nothing
end

@kernel function llg_rhs_stt_cayley_kernel!(dw_dt, @Const(m), @Const(h), @Const(h_stt),
                                            @Const(omega), @Const(pins), alpha::T, beta::T,
                                            gamma::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    if pins[I]
        @inbounds dw_dt[j] = 0
        @inbounds dw_dt[j + 1] = 0
        @inbounds dw_dt[j + 2] = 0
    else
        a = gamma / (1 + alpha * alpha)
        b = alpha * a
        u::T = 1.0 / (1 + alpha * alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz

        @inbounds mht = mx * h_stt[j] + my * h_stt[j + 1] + mz * h_stt[j + 2]
        @inbounds ht1 = mm * h_stt[j] - mht * mx
        @inbounds ht2 = mm * h_stt[j + 1] - mht * my
        @inbounds ht3 = mm * h_stt[j + 2] - mht * mz

        f1 = a * h1 + u * (beta - alpha) * ht1
        f2 = a * h2 + u * (beta - alpha) * ht2
        f3 = a * h3 + u * (beta - alpha) * ht3

        hpt1 = b * h1 + u * (1 + alpha * beta) * ht1
        hpt2 = b * h2 + u * (1 + alpha * beta) * ht2
        hpt3 = b * h3 + u * (1 + alpha * beta) * ht3

        f1 += cross_x(mx, my, mz, hpt1, hpt2, hpt3)
        f2 += cross_y(mx, my, mz, hpt1, hpt2, hpt3)
        f3 += cross_z(mx, my, mz, hpt1, hpt2, hpt3)

        @inbounds wx = omega[j]
        @inbounds wy = omega[j + 1]
        @inbounds wz = omega[j + 2]

        wf = wx * f1 + wy * f2 + wz * f3
        @inbounds dw_dt[j] = f1 - 0.5 * cross_x(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wx
        @inbounds dw_dt[j + 1] = f2 - 0.5 * cross_y(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wy
        @inbounds dw_dt[j + 2] = f3 - 0.5 * cross_z(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wz
    end
end

function llg_rhs_stt_cayley(dw_dt::AbstractArray{T,1}, m::AbstractArray{T,1},
                            h::AbstractArray{T,1}, h_stt::AbstractArray{T,1},
                            omega::AbstractArray{T,1}, alpha::T, beta::T, gamma::T,
                            N::Int64) where {T<:AbstractFloat}
    kernel! = llg_rhs_stt_cayley_kernel!(default_backend[], groupsize[])
    kernel!(dw_dt, m, h, h_stt, omega, alpha, beta, gamma; ndrange=N)
    KernelAbstractions.synchronize(default_backend[])

    return nothing
end

function llg_stt_cayley_call_back(sim::AbstractSim, dw_dt::AbstractArray{T,1}, t::Float64,
                                  omega::AbstractArray{T,1}) where {T<:AbstractFloat}
    driver = sim.driver
    mesh = sim.mesh
    omega_to_spin(omega, sim.prespin, sim.spin, sim.n_total)
    effective_field(sim, sim.spin, t)

    ut = T(driver.ufun(t))
    compute_field_stt(driver.h_stt, sim.spin, driver.ux, driver.uy, driver.uz, mesh.ngbs,
                      ut, T(mesh.dx), T(mesh.dy), T(mesh.dz), sim.n_total)

    llg_rhs_stt_cayley(dw_dt, sim.spin, driver.field, driver.h_stt, omega, sim.pins,
                       driver.alpha, driver.beta, driver.gamma, sim.n_total)

    return nothing
end

@kernel function llg_rhs_stt_cpp_cayley_kernel!(dw_dt, @Const(m), @Const(h), @Const(h_stt),
                                                @Const(omega), @Const(aj), @Const(pins),
                                                px::T, py::T, pz::T, alpha::T, beta::T,
                                                bj::T, gamma::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    if pins[I]
        @inbounds dw_dt[j] = 0
        @inbounds dw_dt[j + 1] = 0
        @inbounds dw_dt[j + 2] = 0
    else
        a = gamma / (1 + alpha * alpha)
        b = alpha * a
        u::T = 1.0 / (1 + alpha * alpha)

        @inbounds aj_ = aj[index]

        @inbounds mx = m[j]
        @inbounds my = m[j + 1]
        @inbounds mz = m[j + 2]
        mm = mx * mx + my * my + mz * mz
        @inbounds mh = mx * h[j] + my * h[j + 1] + mz * h[j + 2]
        @inbounds h1 = mm * h[j] - mh * mx
        @inbounds h2 = mm * h[j + 1] - mh * my
        @inbounds h3 = mm * h[j + 2] - mh * mz

        @inbounds mht = mx * h_stt[j] + my * h_stt[j + 1] + mz * h_stt[j + 2]
        @inbounds ht1 = mm * h_stt[j] - mht * mx
        @inbounds ht2 = mm * h_stt[j + 1] - mht * my
        @inbounds ht3 = mm * h_stt[j + 2] - mht * mz

        mpt = mx * px + my * py + mz * pz
        @inbounds hp1 = mm * px - mpt * mx
        @inbounds hp2 = mm * py - mpt * my
        @inbounds hp3 = mm * pz - mpt * mz

        f1 = a * h1 + u * (beta - alpha) * ht1 + u * (bj - alpha * aj_) * hp1
        f2 = a * h2 + u * (beta - alpha) * ht2 + u * (bj - alpha * aj_) * hp2
        f3 = a * h3 + u * (beta - alpha) * ht3 + u * (bj - alpha * aj_) * hp3

        hpt1 = b * h1 + u * (1 + alpha * beta) * ht1 + u * (aj_ + alpha * bj) * hp1
        hpt2 = b * h2 + u * (1 + alpha * beta) * ht2 + u * (aj_ + alpha * bj) * hp2
        hpt3 = b * h3 + u * (1 + alpha * beta) * ht3 + u * (aj_ + alpha * bj) * hp3

        f1 += cross_x(mx, my, mz, hpt1, hpt2, hpt3)
        f2 += cross_y(mx, my, mz, hpt1, hpt2, hpt3)
        f3 += cross_z(mx, my, mz, hpt1, hpt2, hpt3)

        @inbounds wx = omega[j]
        @inbounds wy = omega[j + 1]
        @inbounds wz = omega[j + 2]

        wf = wx * f1 + wy * f2 + wz * f3
        @inbounds dw_dt[j] = f1 - 0.5 * cross_x(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wx
        @inbounds dw_dt[j + 1] = f2 - 0.5 * cross_y(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wy
        @inbounds dw_dt[j + 2] = f3 - 0.5 * cross_z(wx, wy, wz, f1, f2, f3) + 0.25 * wf * wz
    end
end

function llg_rhs_stt_cpp_cayley(dw_dt::AbstractArray{T,1}, m::AbstractArray{T,1},
                                h::AbstractArray{T,1}, h_stt::AbstractArray{T,1},
                                omega::AbstractArray{T,1}, aj::AbstractArray{T,1},
                                pins::AbstractArray{Bool,1}, p::Tuple, alpha::T, beta::T,
                                bj::T, gamma::T, N::Int64) where {T<:AbstractFloat}
    px, py, pz = T(p[1]), T(p[2]), T(p[3])
    kernel! = llg_rhs_stt_cpp_cayley_kernel!(default_backend[], groupsize[])
    kernel!(dw_dt, m, h, h_stt, omega, aj, pins, px, py, pz, alpha, beta, bj, gamma;
            ndrange=N)
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end

function llg_stt_cpp_cayley_call_back(sim::AbstractSim, dw_dt::AbstractArray{T,1},
                                      t::Float64,
                                      omega::AbstractArray{T,1}) where {T<:AbstractFloat}
    driver = sim.driver
    mesh = sim.mesh
    omega_to_spin(omega, sim.prespin, sim.spin, sim.n_total)
    effective_field(sim, sim.spin, t)
    ut = T(driver.ufun(t))
    compute_field_stt(driver.h_stt, sim.spin, driver.ux, driver.uy, driver.uz, mesh.ngbs,
                      ut, T(mesh.dx), T(mesh.dy), T(mesh.dz), sim.n_total)

    llg_rhs_stt_cpp_cayley(dw_dt, sim.spin, driver.field, driver.h_stt, omega, driver.aj,
                           sim.pins, driver.p, driver.alpha, driver.beta, driver.bj,
                           driver.gamma, sim.n_total)

    return nothing
end
