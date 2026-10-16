#implement the classical Runge-Kutta method for testing purpose

mutable struct RungeKutta{T<:AbstractFloat} <: Integrator
    t::Float64
    step::Float64
    nsteps::Int64
    k1::AbstractArray{T,1}
    k2::AbstractArray{T,1}
    k3::AbstractArray{T,1}
    k4::AbstractArray{T,1}
    y_temp::AbstractArray{T,1}
    dmdt_buf::AbstractArray{T,1}
    max_dmdt::Float64
    rhs_fun::Function
end

function RungeKutta(::Type{T}, n_total::Int64, rhs_fun, step::Float64) where {T<:AbstractFloat}
    k1 = create_zeros(T, 3 * n_total)
    k2 = create_zeros(T, 3 * n_total)
    k3 = create_zeros(T, 3 * n_total)
    k4 = create_zeros(T, 3 * n_total)
    y_temp = create_zeros(T, 3 * n_total)
    dmdt_buf = create_zeros(T, n_total)

    return RungeKutta(0.0, step, 0, k1, k2, k3, k4, y_temp, dmdt_buf, 0.0, rhs_fun)
end

function advance_step(sim::AbstractSim, integrator::RungeKutta)
    h = integrator.step
    k1 = integrator.k1
    k2 = integrator.k2
    k3 = integrator.k3
    k4 = integrator.k4
    y_temp = integrator.y_temp

    #compute k1
    integrator.rhs_fun(sim, k1, sim.spin, integrator.t)

    #compute k2
    y_temp .= sim.spin .+ 0.5 .* h .* k1
    integrator.rhs_fun(sim, k2, y_temp, integrator.t + 0.5 * h)

    #compute k3
    y_temp .= sim.spin .+ 0.5 .* h .* k2
    integrator.rhs_fun(sim, k3, y_temp, integrator.t + 0.5 * h)

    #compute k4
    y_temp .= sim.spin .+ h .* k3
    integrator.rhs_fun(sim, k4, y_temp, integrator.t + h)

    sim.spin .+= (1.0 / 6 * h) .* (k1 .+ 2 .* k2 + 2 .* k3 .+ k4)

    normalise(sim.spin, sim.n_total)
    integrator.nsteps += 1
    integrator.t = integrator.nsteps * h

    # Convergence monitor: k4 is the c = 1 stage, i.e. the RHS at the endpoint
    # predictor — an O(h^2) approximation of the new-state dm/dt (pinned spins
    # are zeroed by the rhs kernel).
    norm! = dmdt_norm_kernel!(get_backend(integrator.dmdt_buf), groupsize[])
    norm!(integrator.dmdt_buf, k4, sim.n_total; ndrange=sim.n_total)
    integrator.max_dmdt = maximum(integrator.dmdt_buf)
    return nothing
end
