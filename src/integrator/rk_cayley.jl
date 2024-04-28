#implement the classical Runge-Kutta method for testing purpose

mutable struct RungeKuttaCayley{T<:AbstractFloat} <: Integrator
    t::Float64
    dt::Float64
    nsteps::Int64
    omega::AbstractArray{T,1}
    dw_dt::AbstractArray{T,1}
    k1::AbstractArray{T,1}
    k2::AbstractArray{T,1}
    k3::AbstractArray{T,1}
    k4::AbstractArray{T,1}
    rhs_fun::Function
end

function RungeKuttaCayley(n_total::Int64, rhs_fun, step::Float64)

    omega = create_zeros(3 * n_total)
    dw_dt = create_zeros(3 * n_total)

    k1 = create_zeros(3 * n_total)
    k2 = create_zeros(3 * n_total)
    k3 = create_zeros(3 * n_total)
    k4 = create_zeros(3 * n_total)

    return RungeKuttaCayley(0.0, step, 0, omega, dw_dt, k1, k2, k3, k4, rhs_fun)
end

function advance_step(sim::AbstractSim, integrator::RungeKuttaCayley)
    h = integrator.dt
    k1 = integrator.k1
    k2 = integrator.k2
    k3 = integrator.k3
    k4 = integrator.k4

    sim.prespin .= sim.spin

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

    omega_to_spin(integrator.omega, sim.prespin, sim.spin, sim.n_total)

    integrator.nsteps += 1
    integrator.t = integrator.nsteps * h
    return nothing
end
