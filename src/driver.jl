mutable struct EmptyDriver <: Driver end

mutable struct EnergyMinimization{T<:AbstractFloat} <: Driver
    gk::AbstractArray{T,1}
    ss::AbstractArray{T,1}
    sf::AbstractArray{T,1}
    ff::AbstractArray{T,1}
    tau::T
    max_tau::T
    min_tau::T
    steps::Int64
end

mutable struct LLG{T<:AbstractFloat} <: Driver
    precession::Bool
    alpha::T
    gamma::T
    ode::Integrator
    tol::Float64
end

mutable struct LLG_CPP{T<:AbstractFloat} <: Driver
    alpha::T
    beta::T
    gamma::T
    aj::AbstractArray{T,1}
    bj::T
    ode::Integrator
    tol::Float64
    p::Tuple{Real,Real,Real}
end

mutable struct LLG_STT{T<:AbstractFloat} <: Driver
    alpha::T
    beta::T
    gamma::T
    ode::Integrator
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
    ode::Integrator
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
    supported_drivers = ["None", "SD", "LLG", "LLG_STT", "LLG_CPP"]
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
        return EnergyMinimization(gk, ss, sf, ff, T(0.0), T(max_tau), T(1e-10), 0)
    end

    supported_integrators = ["Heun", "RungeKutta", "DormandPrince", "DormandPrinceCayley"]
    if !(integrator in supported_integrators)
        error("Supported integrators for GPU: ", join(supported_integrators, " "))
    end

    if driver == "LLG"
        tol = 1e-6
        if integrator == "Heun"
            dopri5 = ModifiedEuler(n_total, llg_call_back, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKutta(n_total, llg_call_back, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrince(n_total, llg_call_back, tol)
        else
            dopri5 = DormandPrinceCayley(n_total, llg_cayley_call_back, tol)
        end
        return LLG(true, T(0.1), T(2.21e5), dopri5, tol)
    elseif driver == "LLG_STT"
        tol = 1e-6
        ux = create_zeros(n_total)
        uy = create_zeros(n_total)
        uz = create_zeros(n_total)
        hstt = KernelAbstractions.zeros(default_backend[], T, 3 * n_total)
        fun = t::Float64 -> 1.0
        if integrator == "Heun"
            dopri5 = ModifiedEuler(n_total, llg_stt_call_back, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKutta(n_total, llg_stt_call_back, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrince(n_total, llg_stt_call_back, tol)
        else
            dopri5 = DormandPrinceCayley(n_total, llg_stt_cayley_call_back, tol)
        end
        return LLG_STT(T(0.5), T(0), T(2.21e5), dopri5, tol, ux, uy, uz, hstt, fun)
    elseif driver == "LLG_CPP"
        tol = 1e-6
        aj = create_zeros(n_total)
        fun = t::Float64 -> 1.0
        p = (0, 0, 1)
        if integrator == "Heun"
            dopri5 = ModifiedEuler(n_total, llg_cpp_call_back, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKutta(n_total, llg_cpp_call_back, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrince(n_total, llg_cpp_call_back, tol)
        else
            error(@sprintf("Unsupported combination driver %s and %s.", driver, integrator))
        end
        return LLG_CPP(T(0.5), T(0), T(2.21e5), aj, T(0), dopri5, tol, p)
    end
    return nothing
end
