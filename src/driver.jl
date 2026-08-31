mutable struct EmptyDriver <: Driver end

mutable struct SD{T<:AbstractFloat} <: Driver
    gk::AbstractArray{T,1}
    ss::AbstractArray{T,1}
    sf::AbstractArray{T,1}
    ff::AbstractArray{T,1}
    tau::T
    max_tau::T
    min_tau::T
    steps::Int64
end

# Landau-Lifshitz-Gilbert driver. The damping `alpha` is an n_total-length array,
# usually an O(1) `Fill` for uniform damping: "spatial or not" is a
# property of the value, not the type, so LLG serves uniform and spatial damping
# alike (the former SpatialLLG driver was merged into it). Use `set_alpha` to
# change it.
mutable struct LLG{T<:AbstractFloat} <: Driver
    precession::Bool
    alpha::AbstractArray{T, 1}
    gamma::T
    integrator::Integrator
    tol::Float64
end

# `sim.driver.alpha = <scalar>` stays supported: the scalar is converted in place to a
# same-length `Fill` (O(1), the length comes from the current array, no mesh needed).
# Arrays pass through with the eltype matched to the driver's precision; functions need
# the mesh for spatial sampling and must go through `set_alpha(sim, ...)` instead.
function Base.setproperty!(driver::LLG, name::Symbol, x)
    if name === :alpha
        T = eltype(driver.alpha)
        n = length(driver.alpha)
        if x isa Function
            throw(ArgumentError("`sim.driver.alpha = <function>` cannot be spatially " *
                                "sampled without the mesh; use `set_alpha(sim, ...)` instead."))
        elseif x isa Number
            x = Fill(T(x), n)
        elseif x isa AbstractArray
            length(x) == n ||
                throw(ArgumentError("`alpha` must have length $n, got $(length(x))."))
            eltype(x) === T || (x = T.(x))
        end
    end
    # default conversion semantics for the other (typed scalar) fields
    return invoke(setproperty!, Tuple{Any,Symbol,Any}, driver, name, x)
end

# Inertial LLG driver (keeps a scalar alpha; not part of the array unification)
mutable struct InertialLLG{T<:AbstractFloat} <: Driver
    alpha::T
    gamma::T
    tau::T  # angular momentum relaxation time
    integrator::Integrator
    tol::Float64
end

# "SpatialLLG" was merged into "LLG" when alpha became a (possibly Fill) array;
# the old name stays as a deprecation alias (remove in v0.7, INTERFACE_DESIGN §10).
function _normalize_driver_name(driver::String)
    if driver == "SpatialLLG"
        Base.depwarn("Driver \"SpatialLLG\" is deprecated; use \"LLG\" together " *
                     "with `set_alpha` (which accepts uniform and spatial alpha).",
                     :create_driver)
        return "LLG"
    end
    return driver
end


function create_driver(driver::String, integrator::String, n_total::Int64)
    if driver in ("LLG_STT", "LLG_CPP", "LLG_STT_CPP")
        replacement = driver == "LLG_CPP" ? "`add_sot(sim; ...)`" :
                      "`add_stt(sim; model=..., ...)`"
        error("The `$driver` driver was removed in v0.6.0. Use the `LLG` driver " *
              "together with $replacement instead.")
    end
    driver = _normalize_driver_name(driver)
    supported_drivers = ["None", "SD", "LLG", "InertialLLG"]
    if !(driver in supported_drivers)
        error("Supported drivers: ", join(supported_drivers, " "))
    end

    T = Float[]

    if driver == "None"
        return EmptyDriver()
    end

    if driver == "SD"
        gk = create_zeros(3 * n_total)
        ss = create_zeros(n_total)
        sf = create_zeros(n_total)
        ff = create_zeros(n_total)
        max_tau = 1.0
        return SD(gk, ss, sf, ff, T(0.0), T(max_tau), T(1e-10), 0)
    end

    supported_integrators = ["Heun", "RungeKutta", "RungeKuttaCayley", "DormandPrince",
                             "DormandPrinceCayley", "DOPRI54", "BS23", "CashKarp54",
                             "Fehlberg54", "RKF54", "GPSM"]
    if !(integrator in supported_integrators)
        error("Supported integrators for GPU: ", join(supported_integrators, " "))
    end

    if driver == "LLG"
        call_back_fun = contains(integrator, "Cayley") ? llg_cayley_call_back :
                        llg_call_back
    elseif driver == "InertialLLG"
        call_back_fun = inertial_llg_call_back
        n_total = 2*n_total
    end

    tol = 1e-6
    if integrator == "Heun"
        dopri5 = ModifiedEuler(n_total, call_back_fun, 1e-14)
    elseif integrator == "RungeKutta"
        dopri5 = RungeKutta(n_total, call_back_fun, 5e-14)
    elseif integrator == "RungeKuttaCayley"
        dopri5 = RungeKuttaCayley(n_total, call_back_fun, 5e-14)
    elseif integrator == "DormandPrince" || integrator == "DOPRI54"
        dopri5 = DormandPrince(n_total, call_back_fun, tol)
    elseif integrator == "BS23"
        dopri5 = BogackiShampine23(n_total, call_back_fun, tol)
    elseif integrator == "CashKarp54"
        dopri5 = CashKarp54(n_total, call_back_fun, tol)
    elseif integrator == "Fehlberg54" || integrator == "RKF54"
        dopri5 = Fehlberg54(n_total, call_back_fun, tol)
    elseif integrator == "DormandPrinceCayley"
        dopri5 = DormandPrinceCayley(n_total, call_back_fun, tol)
    elseif integrator == "GPSM"
        dopri5 = GPSM(n_total, 1e-13)
    end

    if driver == "LLG"
        # the default uniform alpha is an O(1) Fill: zero allocation
        return LLG(true, Fill(T(0.1), n_total), T(2.21e5), dopri5, tol)
    elseif driver == "InertialLLG"
        return InertialLLG(T(0.01), T(2.21e5), T(10e-12), dopri5, tol)
    end

    return nothing
end
