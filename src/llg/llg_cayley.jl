
# `alphas` is read per spin so uniform (Fill) and spatial (dense) damping share
# this kernel, mirroring llg_rhs_kernel!.
@kernel function llg_rhs_cayley_kernel!(dw_dt, @Const(m), @Const(h), @Const(omega),
                                        pins, alphas, gamma::T,
                                        precession::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds pin::Bool = pins[I]
    @inbounds alpha::T = alphas[I]

    if pin
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
        @inbounds dw_dt[j] = f1 - T(0.5) * cross_x(wx, wy, wz, f1, f2, f3) + T(0.25) * wf * wx
        @inbounds dw_dt[j + 1] = f2 - T(0.5) * cross_y(wx, wy, wz, f1, f2, f3) + T(0.25) * wf * wy
        @inbounds dw_dt[j + 2] = f3 - T(0.5) * cross_z(wx, wy, wz, f1, f2, f3) + T(0.25) * wf * wz
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

    kernel! = llg_rhs_cayley_kernel!(get_backend(sim.spin), groupsize[])
    kernel!(dw_dt, sim.spin, sim.field, omega, sim.pins, driver.alpha, driver.gamma,
            driver.precession; ndrange=N)

    return nothing
end

