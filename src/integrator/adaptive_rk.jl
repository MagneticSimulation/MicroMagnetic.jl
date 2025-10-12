# Description: Implementation of adaptive Runge-Kutta integrators for ODEs

abstract type AbstractRKMethod end

struct DOPRI54Method <: AbstractRKMethod
    name::String
    stages::Int
    order_main::Int
    order_embedded::Int
    FSAL::Bool
end

struct BS23Method <: AbstractRKMethod
    name::String
    stages::Int
    order_main::Int
    order_embedded::Int
    FSAL::Bool
end

struct CashKarp54Method <: AbstractRKMethod
    name::String
    stages::Int
    order_main::Int
    order_embedded::Int
    FSAL::Bool
end

struct Fehlberg54Method <: AbstractRKMethod
    name::String
    stages::Int
    order_main::Int
    order_embedded::Int
    FSAL::Bool
end

mutable struct AdaptiveRK{T<:AbstractFloat,A<:AbstractArray{T}} <: Integrator
    method::AbstractRKMethod
    tol::Float64
    t::Float64
    step::Float64
    step_next::Float64
    facmax::Float64
    facmin::Float64
    safety::Float64
    nsteps::Int64
    nfevals::Int64
    y_current::A
    y_next::A
    ks::Vector{A}
    rhs_fun::Function
    succeed::Bool
    step_function::Function  # Method-specific step function
    y_temp::A # Pre-allocated temporary storage for RK stages
end

function AdaptiveRK(n_total::Int64, method::M, rhs_fun,
                    tol::Float64) where {M<:AbstractRKMethod}
    stages = method.stages

    # Allocate storage
    y_current = create_zeros(3 * n_total)
    y_next = create_zeros(3 * n_total)
    ks = [create_zeros(3 * n_total) for _ in 1:stages]
    y_temp = create_zeros(3 * n_total)  # Pre-allocate temporary storage

    facmax = 5.0
    facmin = 0.2
    safety = 0.824

    # Select method-specific step function
    step_function = if method isa DOPRI54Method
        rk_step_dopri54!
    elseif method isa BS23Method
        rk_step_bs23!
    elseif method isa CashKarp54Method
        rk_step_cashkarp54!
    elseif method isa Fehlberg54Method
        rk_step_fehlberg54!
    else
        error("Unknown RK method")
    end

    return AdaptiveRK(method, tol, 0.0, 0.0, 0.0, facmax, facmin, safety, 0, 0, y_current,
                      y_next, ks, rhs_fun, false, step_function, y_temp)
end

# Generic step function (now delegates to method-specific implementation)
function rk_step_inner(sim::AbstractSim, step::Float64, t::Float64, integrator::AdaptiveRK)
    return integrator.step_function(sim, step, t, integrator)
end

function advance_step(sim::AbstractSim, integrator::AdaptiveRK)
    max_nsteps = 50

    t = integrator.t
    order = integrator.method.order_embedded

    # Initialize step size if needed
    if integrator.step_next <= 0
        integrator.step_next = compute_initial_step(sim, integrator, 1e-12)
    end

    integrator.step = integrator.step_next

    nstep = 1
    while true
        # Perform the RK step
        max_error = rk_step_inner(sim, integrator.step, t, integrator) / integrator.tol

        # Handle NaN errors
        if isnan(max_error)
            integrator.step_next = 1e-14
            @warn "NaN error detected, reducing step size"
            return false
        end

        nstep += 1
        if nstep > max_nsteps
            @warn "Too many inner steps ($nstep), system may not converge!"
            return false
        end

        integrator.succeed = (max_error <= 1.0)

        if integrator.succeed
            # Update state
            integrator.y_current .= integrator.y_next

            # Update statistics
            integrator.nsteps += 1
            integrator.t += integrator.step

            # Calculate next step size
            factor = integrator.safety * (1.0 / max_error)^(1/(order+1))
            integrator.step_next = integrator.step *
                                   min(integrator.facmax, max(integrator.facmin, factor))
            break
        else
            # Reduce step size and retry
            factor = integrator.safety * (1.0 / max_error)^(1/order)
            integrator.step = integrator.step *
                              min(integrator.facmax, max(integrator.facmin, factor))
        end
    end

    return true
end

# Simple initial step size estimation
function compute_initial_step(sim::AbstractSim, integrator::AdaptiveRK, dt)
    abs_step = dt
    abs_step_tmp = dt
    integrator.step = 1e-15
    integrator.rhs_fun(sim, integrator.y_temp, integrator.y_current, integrator.t)
    order = integrator.method.order_main
    r_step = maximum(abs.(integrator.y_temp)) /
             (integrator.safety * integrator.tol^(1 / order))
    integrator.nfevals += 1
    #TODO: find a better init step?
    if abs_step * r_step > 0.001
        abs_step_tmp = 0.001 / r_step
    end
    return min(abs_step, abs_step_tmp)
end

# Additional utility functions for time stepping
function integrate_to_time(sim::AbstractSim, integrator::AdaptiveRK, t_final::Float64)

    # Integrate from current time to t_final
    while integrator.t < t_final
        # Adjust step to not overshoot final time
        if integrator.step_next + integrator.t > t_final
            integrator.step_next = t_final - integrator.t
        end

        success = advance_step(sim, integrator)

        if !success
            @error "Integration failed at t = $(integrator.t)"
            break
        end
    end

    return integrator.succeed
end

function get_statistics(integrator::AdaptiveRK)
    # Return integration statistics
    return (time=integrator.t, steps=integrator.nsteps, nfevals=integrator.nfevals,
            current_step=integrator.step, success=integrator.succeed)
end

function reset_integrator!(integrator::AdaptiveRK)
    # Reset integrator to initial state
    integrator.t = 0.0
    integrator.step = 0.0
    integrator.step_next = 0.0
    integrator.nsteps = 0
    integrator.nfevals = 0
    integrator.succeed = false

    # Reset all arrays to zero
    integrator.y_current .= 0.0
    integrator.y_next .= 0.0
    integrator.y_temp .= 0.0
    for k in integrator.ks
        k .= 0.0
    end
end

# Initialization functions
function initialize!(integrator::AdaptiveRK, y0::AbstractArray)
    integrator.y_current .= y0
    integrator.t = 0.0
    integrator.step = 0.0
    integrator.step_next = 0.0
    integrator.nsteps = 0
    integrator.nfevals = 0
    integrator.succeed = false

    # Reset all k vectors and temporary storage to zero
    for k in integrator.ks
        k .= 0.0
    end
    return integrator.y_temp .= 0.0  # Also reset temporary storage
end

function set_initial_condition!(integrator::AdaptiveRK, y0::AbstractArray)
    return integrator.y_current .= y0
end

function get_current_state(integrator::AdaptiveRK)
    return integrator.y_current
end

function get_current_time(integrator::AdaptiveRK)
    return integrator.t
end

# Constructors for specific methods
function DormandPrince(n_total::Int64, rhs_fun, tol::Float64)
    method = DOPRI54Method("DOPRI54", 7, 5, 4, false)
    return AdaptiveRK(n_total, method, rhs_fun, tol)
end

function BogackiShampine23(n_total::Int64, rhs_fun, tol::Float64)
    method = BS23Method("BS23", 4, 3, 2, true)
    return AdaptiveRK(n_total, method, rhs_fun, tol)
end

function CashKarp54(n_total::Int64, rhs_fun, tol::Float64)
    method = CashKarp54Method("CashKarp54", 6, 5, 4, false)
    return AdaptiveRK(n_total, method, rhs_fun, tol)
end

function Fehlberg54(n_total::Int64, rhs_fun, tol::Float64)
    method = Fehlberg54Method("Fehlberg54", 6, 5, 4, false)
    return AdaptiveRK(n_total, method, rhs_fun, tol)
end

function run_step(sim::AbstractSim, driver::LLG)
    return advance_step(sim, driver.integrator)
end

function run_step(sim::AbstractSim, driver::LLG_STT)
    return advance_step(sim, driver.integrator)
end
