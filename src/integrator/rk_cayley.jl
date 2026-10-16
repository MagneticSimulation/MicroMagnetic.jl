#implement the classical Runge-Kutta method for testing purpose

mutable struct RungeKuttaCayley{T<:AbstractFloat} <: Integrator
    t::Float64
    step::Float64
    nsteps::Int64
    omega::AbstractArray{T,1}
    dw_dt::AbstractArray{T,1}
    k1::AbstractArray{T,1}
    k2::AbstractArray{T,1}
    k3::AbstractArray{T,1}
    k4::AbstractArray{T,1}
    m_pre::AbstractArray{T,1}       # spin at the start of the current step
    dm_dt::AbstractArray{T,1}       # scratch for the convergence monitor
    dmdt_buf::AbstractArray{T,1}    # per-cell |dm/dt| scratch
    max_dmdt::Float64               # max |dm/dt| of the last step (rad/s)
    rhs_fun::Function
end

function RungeKuttaCayley(::Type{T}, n_total::Int64, rhs_fun, step::Float64) where {T<:AbstractFloat}

    omega = create_zeros(T, 3 * n_total)
    dw_dt = create_zeros(T, 3 * n_total)

    k1 = create_zeros(T, 3 * n_total)
    k2 = create_zeros(T, 3 * n_total)
    k3 = create_zeros(T, 3 * n_total)
    k4 = create_zeros(T, 3 * n_total)
    m_pre = create_zeros(T, 3 * n_total)
    dm_dt = create_zeros(T, 3 * n_total)
    dmdt_buf = create_zeros(T, n_total)

    return RungeKuttaCayley(0.0, step, 0, omega, dw_dt, k1, k2, k3, k4, m_pre, dm_dt,
                            dmdt_buf, 0.0, rhs_fun)
end

function advance_step(sim::AbstractSim, integrator::RungeKuttaCayley)
    h = integrator.step
    k1 = integrator.k1
    k2 = integrator.k2
    k3 = integrator.k3
    k4 = integrator.k4

    integrator.m_pre .= sim.spin

    y_next = integrator.omega
    fill!(y_next, 0) # we always have y=0

    #compute k1
    integrator.rhs_fun(sim, k1, integrator.t, y_next)

    #compute k2
    y_next .= 0.5 .* h .* k1
    integrator.rhs_fun(sim, k2, integrator.t + 0.5 * h, y_next)

    #compute k3
    y_next .= 0.5 .* h .* k2
    integrator.rhs_fun(sim, k3, integrator.t + 0.5 * h, y_next)

    #compute k4
    y_next  .=  h .* k3
    integrator.rhs_fun(sim, k4, integrator.t + h, y_next)

    y_next .= (1.0 / 6 * h) .* (k1 .+ 2 .* k2 + 2 .* k3 .+ k4)

    omega_to_spin(integrator.omega, integrator.m_pre, sim.spin, sim.n_total)

    integrator.nsteps += 1
    integrator.t = integrator.nsteps * h

    # Convergence monitor: full LLG RHS at the new state. The field left by the
    # last stage call corresponds to the c = 1 predictor state (an O(h^2)
    # offset), which is plenty for a stopping criterion. Pinned spins are
    # zeroed by the kernel.
    driver = sim.driver
    T = eltype(integrator.dm_dt)
    kernel! = llg_rhs_kernel!(get_backend(sim.spin), groupsize[])
    kernel!(integrator.dm_dt, sim.spin, sim.field, sim.pins, driver.alpha,
            T(driver.gamma), driver.precession; ndrange=sim.n_total)
    norm! = dmdt_norm_kernel!(get_backend(integrator.dmdt_buf), groupsize[])
    norm!(integrator.dmdt_buf, integrator.dm_dt, sim.n_total; ndrange=sim.n_total)
    integrator.max_dmdt = maximum(integrator.dmdt_buf)
    return nothing
end
