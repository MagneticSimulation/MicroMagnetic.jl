#implement the  Runge-Kutta method with adaptive stepsize using Dormand–Prince pair
#https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Dormand%E2%80%93Prince

#To use this integrator, one needs to provide:
#   - A simulation instance with field "spin"
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
    y_next::AbstractArray{T}
    k1::AbstractArray{T}
    k2::AbstractArray{T}
    k3::AbstractArray{T}
    k4::AbstractArray{T}
    k5::AbstractArray{T}
    k6::AbstractArray{T}
    k7::AbstractArray{T}
    rhs_fun::Function
    succeed::Bool
    FSAL::Bool
end

function DormandPrince(n_total::Int64, rhs_fun, tol::Float64)
    T = Float[]
    y_next = create_zeros(3 * n_total)
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
    return DormandPrince(tol, 0.0, 0.0, 0.0, facmax, facmin, safety, 0, 0, y_next, k1, k2,
                         k3, k4, k5, k6, k7, rhs_fun, false, true)
end

function dopri5_step_inner(sim::AbstractSim, step::Float64, t::Float64)
    a = (1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0)
    b = (1 / 5, 3 / 40, 9 / 40, 44 / 45, -56 / 15, 32 / 9)
    c = (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729)
    d = (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656)
    v = (35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84)
    w = (71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)
    I = sim.driver.integrator

    k1, k2, k3, k4, k5, k6, k7 = I.k1, I.k2, I.k3, I.k4, I.k5, I.k6, I.k7
    y_next = I.y_next
    y_current = sim.spin

    if I.FSAL && I.succeed && t >= I.t
        k1 .= k7 # copy k7 directly to k1
    else
        I.rhs_fun(sim, k1, y_current, t) #compute k1
        I.nfevals += 1
    end

    # compute k2
    vector_add2(y_next, y_current, k1, b[1]*step)
    I.rhs_fun(sim, k2, y_next, t + a[1] * step) 

    # compute k3
    vector_add3(y_next, y_current, k1, k2, b[2]*step, b[3]*step)
    I.rhs_fun(sim, k3, y_next, t + a[2] * step) 

    # compute k4
    vector_add4(y_next, y_current, k1, k2, k3, b[4]*step, b[5]*step, b[6]*step)
    I.rhs_fun(sim, k4, y_next, t + a[3] * step) 

    # #compute k5
    vector_add5(y_next, y_current, k1, k2, k3, k4, c[1]*step, c[2]*step, c[3]*step, c[4]*step)
    I.rhs_fun(sim, k5, y_next, t + a[4] * step) 

    # compute k6
    vector_add6(y_next, y_current, k1, k2, k3, k4, k5, d[1]*step, d[2]*step, d[3]*step, d[4]*step, d[5]*step)
    I.rhs_fun(sim, k6, y_next, t + a[5] * step) 

    # compute k7. note v[2] = 0, so we still use vector_add6
    vector_add6(y_next, y_current, k1, k3, k4, k5, k6, v[1]*step, v[3]*step, v[4]*step, v[5]*step, v[6]*step)
    normalise(y_next, sim.n_total) #if we want to copy k7 to k1, we should normalise it here.
    I.rhs_fun(sim, k7, y_next, t + a[6] * step) 

    I.nfevals += 6
    # compute the error, we use k2 to store the error since w[2] = 0
    vector_add6b(k2, k1, k3, k4, k5, k6, k7, w[1], w[3], w[4], w[5], w[6], w[7])

    max_error = maximum(abs.(k2))*step + eps()

    return max_error
end

#this function can be run in both GPU and CPU
function compute_init_step_DP(sim::AbstractSim, dt::Float64)
    abs_step = dt
    abs_step_tmp = dt
    integrator = sim.driver.integrator
    integrator.step = 1e-15
    integrator.rhs_fun(sim, integrator.k2, sim.spin, integrator.t)
    #abs!(integrator.errors)
    r_step = maximum(abs.(integrator.k2)) / (integrator.safety * integrator.tol^0.2)
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

    # TODO: delete sim.prespin since we have y_next
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
            sim.spin .= integrator.y_next
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
