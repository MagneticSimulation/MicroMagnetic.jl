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

# Landau-Lifshitz-Gilbert driver
mutable struct LLG{T<:AbstractFloat} <: Driver
    precession::Bool
    alpha::T
    gamma::T
    integrator::Integrator
    tol::Float64
end

# Inertial LLG driver
mutable struct InertialLLG{T<:AbstractFloat} <: Driver
    alpha::T
    gamma::T
    tau::T  # angular momentum relaxation time
    integrator::Integrator
    tol::Float64
end

mutable struct SpatialLLG{T<:AbstractFloat} <: Driver
    precession::Bool
    alpha::AbstractArray{T, 1}
    gamma::T
    integrator::Integrator
    tol::Float64
end

mutable struct LLG_CPP{T<:AbstractFloat} <: Driver
    alpha::T
    beta::T
    gamma::T
    aj::AbstractArray{T,1}
    bj::T
    integrator::Integrator
    tol::Float64
    p::Tuple{Real,Real,Real}
    ufun::Function
end

mutable struct LLG_STT{T<:AbstractFloat} <: Driver
    alpha::T
    beta::T
    gamma::T
    integrator::Integrator
    tol::Float64
    ux::AbstractArray{T,1}
    uy::AbstractArray{T,1}
    uz::AbstractArray{T,1}
    h_stt::AbstractArray{T,1}
    ufun::Function
end

mutable struct LLG_STT_CPP{T<:AbstractFloat} <: Driver
    alpha::T
    beta::T
    gamma::T
    bj::T
    integrator::Integrator
    tol::Float64
    p::Tuple{Real,Real,Real}
    aj::AbstractArray{T,1}
    ux::AbstractArray{T,1}
    uy::AbstractArray{T,1}
    uz::AbstractArray{T,1}
    h_stt::AbstractArray{T,1}
    ufun::Function
end

function create_driver(driver::String, integrator::String, n_total::Int64)
    supported_drivers = ["None", "SD", "LLG", "LLG_STT", "LLG_CPP", "SpatialLLG", "InertialLLG"]
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
                             "Fehlberg54", "RKF54"]
    if !(integrator in supported_integrators)
        error("Supported integrators for GPU: ", join(supported_integrators, " "))
    end

    if driver == "LLG"
        call_back_fun = contains(integrator, "Cayley") ? llg_cayley_call_back :
                        llg_call_back
    elseif driver == "LLG_STT"
        @warn "The driver 'LLG_STT' is deprecated. Please use 'add_stt' function. The 'LLG_STT' will be removed in v0.5.0"
        call_back_fun = contains(integrator, "Cayley") ? llg_stt_cayley_call_back :
                        llg_stt_call_back
    elseif driver == "LLG_CPP"
        @warn "The driver 'LLG_CPP' is deprecated. Please use 'add_sot' function. The 'LLG_CPP' will be removed in v0.5.0"
        call_back_fun = llg_cpp_call_back
    elseif driver == "SpatialLLG"
        call_back_fun = spatial_llg_call_back
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
    end

    if driver == "LLG"
        return LLG(true, T(0.1), T(2.21e5), dopri5, tol)
    elseif driver == "SpatialLLG"
        alpha = create_zeros(n_total)
        alpha .= 0.1
        return SpatialLLG(true, alpha, T(2.21e5), dopri5, tol)
    elseif driver == "InertialLLG"
        return InertialLLG(T(0.01), T(2.21e5), T(10e-12), dopri5, tol)
    elseif driver == "LLG_STT"
        tol = 1e-6
        ux = create_zeros(n_total)
        uy = create_zeros(n_total)
        uz = create_zeros(n_total)
        hstt = KernelAbstractions.zeros(default_backend[], T, 3 * n_total)
        fun = t::Float64 -> 1.0
        return LLG_STT(T(0.5), T(0), T(2.21e5), dopri5, tol, ux, uy, uz, hstt, fun)
    elseif driver == "LLG_CPP"
        tol = 1e-6
        aj = create_zeros(n_total)
        fun = t::Float64 -> 1.0
        p = (0, 0, 1)
        if integrator == "DormandPrinceCayley"
            error(@sprintf("Unsupported combination driver %s and %s.", driver, integrator))
        end
        return LLG_CPP(T(0.5), T(0), T(2.21e5), aj, T(0), dopri5, tol, p, fun)
    end
    return nothing
end
