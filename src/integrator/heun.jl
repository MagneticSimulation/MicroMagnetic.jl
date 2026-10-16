#implement the Heun (Modified Euler) algorithm,
#https://en.wikipedia.org/wiki/Heun%27s_method

# k1 = f(t_n, y_n)
# k2 = f(t_n+1, y_n + h*k1)
# y_n+1 = y_n + 1/2*h*(k1+k2)

mutable struct ModifiedEuler{T<:AbstractFloat} <: Integrator
   t::Float64
   step::Float64
   nsteps::Int64
   k1::AbstractArray{T,1}
   k2::AbstractArray{T,1}
   y_temp::AbstractArray{T,1}
   dmdt_buf::AbstractArray{T,1}
   max_dmdt::Float64
   rhs_fun::Function
end

function ModifiedEuler(::Type{T}, n_total::Int64, rhs_fun, step::Float64) where {T<:AbstractFloat}
  k1 = create_zeros(T, 3 * n_total)
  k2 = create_zeros(T, 3 * n_total)
  y_temp = create_zeros(T, 3 * n_total)
  dmdt_buf = create_zeros(T, n_total)
  return ModifiedEuler(0.0, step, 0, k1, k2, y_temp, dmdt_buf, 0.0, rhs_fun)
end

function advance_step(sim::AbstractSim, integrator::ModifiedEuler)
    h = integrator.step
    k1 =  integrator.k1
    k2 =  integrator.k2
    y_temp = integrator.y_temp

    #compute k1
    integrator.rhs_fun(sim, k1, sim.spin, integrator.t)

    #compute k2
    vector_add2(y_temp, sim.spin, k1, h)
    integrator.rhs_fun(sim, k2, y_temp, integrator.t+h)

    sim.spin .+= 0.5 .* h .* (k1 .+ k2)

    normalise(sim.spin, sim.n_total)
    integrator.nsteps += 1
    integrator.t = integrator.nsteps*h

    # Convergence monitor: k2 is the c = 1 stage, i.e. the RHS at the endpoint
    # predictor — an O(h^2) approximation of the new-state dm/dt (pinned spins
    # are zeroed by the rhs kernel).
    norm! = dmdt_norm_kernel!(get_backend(integrator.dmdt_buf), groupsize[])
    norm!(integrator.dmdt_buf, k2, sim.n_total; ndrange=sim.n_total)
    integrator.max_dmdt = maximum(integrator.dmdt_buf)

    return nothing

end
