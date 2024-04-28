#implement the  Runge-Kutta method with adaptive stepsize using Dormand–Prince pair
#https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Dormand%E2%80%93Prince

#To use this integrator, one needs to provide:
#   - A simulation instance with field "spin" and "prespin"
#   - A call back function

mutable struct DormandPrince{T<:AbstractFloat} <: Integrator
    tol::Float64
    t::Float64
    step::Float64
    step_next::Float64
    facmax::Float64
    facmin::Float64
    safety::Float64
    nsteps::Int64
    nfevals::Int64
    errors::AbstractArray{T}
    k1::AbstractArray{T}
    k2::AbstractArray{T}
    k3::AbstractArray{T}
    k4::AbstractArray{T}
    k5::AbstractArray{T}
    k6::AbstractArray{T}
    k7::AbstractArray{T}
    rhs_fun::Function
    succeed::Bool
end

function DormandPrince(n_total::Int64, rhs_fun, tol::Float64)
    T = Float[]
    errors = create_zeros(3 * n_total)
    k1 = create_zeros(3 * n_total)
    k2 = create_zeros(3 * n_total)
    k3 = create_zeros(3 * n_total)
    k4 = create_zeros(3 * n_total)
    k5 = create_zeros(3 * n_total)
    k6 = create_zeros(3 * n_total)
    k7 = create_zeros(3 * n_total)
    facmax = 5.0
    facmin = 0.2
    safety = 0.824
    return DormandPrince(tol, 0.0, 0.0, 0.0, facmax, facmin, safety, 0, 0, errors, k1, k2,
                         k3, k4, k5, k6, k7, rhs_fun, false)
end

#this function can be run in both GPU and CPU, however, does not work for GPU+MPI?
function dopri5_step_inner(sim::AbstractSim, step::Float64, t::Float64)
    a = (1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0)
    b = (1 / 5, 3 / 40, 9 / 40, 44 / 45, -56 / 15, 32 / 9)
    c = (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729)
    d = (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656)
    v = (35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84)
    w = (71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)
    ode = sim.driver.integrator
    rhs = ode.rhs_fun

    k1, k2, k3, k4, k5, k6, k7 = ode.k1, ode.k2, ode.k3, ode.k4, ode.k5, ode.k6, ode.k7

    y_next = sim.spin
    y_current = sim.prespin
    y_next .= y_current
    #println("From integrator: y_next", y_next)
    ode.rhs_fun(sim, k1, y_next, t) #compute k1, TODO: copy k7 directly to k1

    y_next .= y_current .+ b[1] .* k1 .* step
    ode.rhs_fun(sim, k2, y_next, t + a[1] * step) #k2

    y_next .= y_current .+ (b[2] .* k1 .+ b[3] .* k2) .* step
    ode.rhs_fun(sim, k3, y_next, t + a[2] * step) #k3

    y_next .= y_current .+ (b[4] .* k1 .+ b[5] .* k2 .+ b[6] .* k3) .* step
    ode.rhs_fun(sim, k4, y_next, t + a[3] * step) #k4

    y_next .= y_current .+ (c[1] .* k1 .+ c[2] .* k2 + c[3] .* k3 .+ c[4] .* k4) .* step
    ode.rhs_fun(sim, k5, y_next, t + a[4] * step) #k5

    y_next .= y_current .+
              (d[1] .* k1 .+ d[2] .* k2 .+ d[3] .* k3 .+ d[4] .* k4 + d[5] .* k5) .* step
    ode.rhs_fun(sim, k6, y_next, t + a[5] * step) #k6

    y_next .= y_current .+
              (v[1] .* k1 .+ v[2] .* k2 .+ v[3] .* k3 .+ v[4] .* k4 .+ v[5] .* k5 .+
               v[6] .* k6) .* step
    normalise(y_next, sim.n_total) #if we want to copy k7 to k1, we should normalise it here.
    ode.rhs_fun(sim, k7, y_next, t + a[6] * step) #k7

    ode.nfevals += 7
    ode.errors .= (w[1] .* k1 + w[2] .* k2 .+ w[3] .* k3 .+ w[4] .* k4 .+ w[5] .* k5 +
                   w[6] .* k6 +
                   w[7] .* k7) .* step

    max_error = maximum(abs.(ode.errors)) + eps() #TODO: check whether CUDA has fixed this eat memory bug
    #abs!(ode.errors)
    #max_error =  maximum(ode.errors) + eps()

    return max_error
end

#this function can be run in both GPU and CPU
function compute_init_step_DP(sim::AbstractSim, dt::Float64)
    abs_step = dt
    abs_step_tmp = dt
    integrator = sim.driver.integrator
    integrator.step = 1e-15
    integrator.rhs_fun(sim, integrator.errors, sim.spin, integrator.t)
    #abs!(integrator.errors)
    r_step = maximum(abs.(integrator.errors)) / (integrator.safety * integrator.tol^0.2)
    integrator.nfevals += 1
    #FIXME: how to obtain a reasonable init step?
    if abs_step * r_step > 0.001
        abs_step_tmp = 0.001 / r_step
    end
    return min(abs_step, abs_step_tmp)
end

function advance_step(sim::AbstractSim, integrator::DormandPrince)
    max_nsteps = 100

    t = integrator.t

    sim.prespin .= sim.spin

    if integrator.step_next <= 0
        integrator.step_next = compute_init_step_DP(sim, 1e-12)
        if integrator.step_next < 1e-15
            #integrator.step_next = 1e-15
        end
    end

    integrator.step = integrator.step_next

    nstep = 1
    while true
        max_error = dopri5_step_inner(sim, integrator.step, t) / integrator.tol
        if isnan(max_error)
            step_next = 1e-14
        end
        nstep += 1
        if nstep > max_nsteps
            @info("Too many inner steps, the system may not converge!", step_next,
                  max_error)
            return false
        end
        integrator.succeed = (max_error <= 1)

        if integrator.succeed
            integrator.nsteps += 1
            #integrator.step = step_next
            integrator.t += integrator.step
            factor = integrator.safety * (1.0 / max_error)^0.2
            integrator.step_next = integrator.step *
                                   min(integrator.facmax, max(integrator.facmin, factor))
            break
        else
            factor = integrator.safety * (1.0 / max_error)^0.25
            integrator.step = integrator.step *
                              min(integrator.facmax, max(integrator.facmin, factor))
        end
    end
    return true
end

function run_step(sim::AbstractSim, driver::LLG)
    return advance_step(sim, driver.integrator)
end

function run_step(sim::AbstractSim, driver::LLG_STT)
    return advance_step(sim, driver.integrator)
end
