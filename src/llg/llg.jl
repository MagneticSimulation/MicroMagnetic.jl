"""
GPU kernel to compute the standard LLG equation,
    dm/dt = - gamma_L * (m x H) - alpha*gamma_L* m x (m x H)
where gamma_L = gamma/(1+alpha^2).

The damping `alphas` is read per spin (`alphas[I]`), so the same kernel serves a
dense spatial-alpha array and an O(1) `Fill` for uniform damping; for an isbits
Fill the read folds into an immediate.
"""
@kernel function llg_rhs_kernel!(dm_dt, @Const(m), @Const(h), pins, alphas,
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
function llg_call_back(sim::AbstractSim, dm_dt::AbstractArray{T,1},
                       spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    driver = sim.driver

    effective_field(sim, spin, t)

    @timeit MicroMagnetic.timer "llg" begin
        kernel! = llg_rhs_kernel!(get_backend(sim.spin), groupsize[])
        kernel!(dm_dt, spin, sim.field, sim.pins, driver.alpha, T(driver.gamma),
                driver.precession; ndrange=N)
    end

    return nothing
end
