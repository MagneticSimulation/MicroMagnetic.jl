#implement the classical Runge-Kutta method for testing purpose

mutable struct RungeKutta{T<:AbstractFloat} <: Integrator
    t::Float64
    step::Float64
    nsteps::Int64
    k1::AbstractArray{T,1}
    k2::AbstractArray{T,1}
    k3::AbstractArray{T,1}
    k4::AbstractArray{T,1}
    rhs_fun::Function
end

function RungeKutta(n_total::Int64, rhs_fun, step::Float64)
    T = single_precision.x ? Float32 : Float64
    k1 = create_zeros(3 * n_total)
    k2 = create_zeros(3 * n_total)
    k3 = create_zeros(3 * n_total)
    k4 = create_zeros(3 * n_total)

    return RungeKutta(0.0, step, 0, k1, k2, k3, k4, rhs_fun)
end

function advance_step(sim::AbstractSim, integrator::RungeKutta)
    h = integrator.step
    k1 = integrator.k1
    k2 = integrator.k2
    k3 = integrator.k3
    k4 = integrator.k4

    #compute k1
    integrator.rhs_fun(sim, k1, sim.spin, integrator.t)

    #compute k2
    sim.prespin .= sim.spin .+ 0.5 .* h .* k1
    integrator.rhs_fun(sim, k2, sim.prespin, integrator.t + 0.5 * h)

    #compute k3
    sim.prespin .= sim.spin .+ 0.5 .* h .* k2
    integrator.rhs_fun(sim, k3, sim.prespin, integrator.t + 0.5 * h)

    #compute k4
    sim.prespin .= sim.spin .+ h .* k3
    integrator.rhs_fun(sim, k4, sim.prespin, integrator.t + h)

    sim.prespin .= sim.spin
    sim.spin .+= (1.0 / 6 * h) .* (k1 .+ 2 .* k2 + 2 .* k3 .+ k4)

    normalise(sim.spin, sim.n_total)
    integrator.nsteps += 1
    return integrator.t = integrator.nsteps * h
end
