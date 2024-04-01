mutable struct EmptyDriverGPU <: DriverGPU
end

mutable struct EnergyMinimizationGPU{T<:AbstractFloat} <: DriverGPU
  gk::CuArray{T, 1}
  ss::CuArray{T, 1}
  sf::CuArray{T, 1}
  ff::CuArray{T, 1}
  field::CuArray{T, 1}
  tau::T
  max_tau::T
  min_tau::T
  steps::Int64
end

mutable struct LLG_GPU{T<:AbstractFloat} <: DriverGPU
  precession::Bool
  alpha::T
  gamma::T
  ode::Integrator
  tol::Float64
  field::CuArray{T, 1}
end

mutable struct LLG_STT_GPU{T<:AbstractFloat} <: DriverGPU
  alpha::T
  beta::T
  gamma::T
  ode::Integrator
  tol::Float64
  ux::CuArray{T, 1}
  uy::CuArray{T, 1}
  uz::CuArray{T, 1}
  h_stt::CuArray{T, 1}
  field::CuArray{T, 1}
  ufun::Function
end

mutable struct LLG_CPP_GPU{T<:AbstractFloat} <: DriverGPU
  alpha::T
  beta::T
  gamma::T
  aj::CuArray{T, 1}
  bj::T
  ode::Integrator
  tol::Float64
  p::Tuple{Real, Real, Real}
  field::CuArray{T, 1}
end

mutable struct LLG_STT_CPP_GPU{T<:AbstractFloat} <: DriverGPU
  alpha::T
  beta::T
  gamma::T
  bj::T
  ode::Integrator
  tol::Float64
  p::Tuple{Real, Real, Real}
  aj::CuArray{T, 1}
  ux::CuArray{T, 1}
  uy::CuArray{T, 1}
  uz::CuArray{T, 1}
  h_stt::CuArray{T, 1}
  field::CuArray{T, 1}
  ufun::Function
end

function create_driver_gpu(driver::String, integrator::String, n_total::Int64)
    supported_drivers = ["None", "SD", "LLG", "LLG_STT", "LLG_CPP", "LLG_STT_CPP"]
    if !(driver in supported_drivers)
        error("Supported drivers for GPU: ", join(supported_drivers, " "))
    end

    Float = _cuda_using_double.x ? Float64 : Float32

    if driver == "None"
        return EmptyDriverGPU()
    end

    if driver == "SD"
        gk = CUDA.zeros(Float,3*n_total)
        field = CUDA.zeros(Float,3*n_total)
        ss = CUDA.zeros(Float, n_total)
        sf = CUDA.zeros(Float, n_total)
        ff = CUDA.zeros(Float, n_total)
        max_tau = 1.0
        return EnergyMinimizationGPU(gk, ss, sf, ff, field, Float(0.0), Float(max_tau), Float(1e-10), 0)
    end

    supported_integrators = ["Heun", "RungeKutta", "DormandPrince", "DormandPrinceCayley"]
    if !(integrator in supported_integrators)
        error("Supported integrators for GPU: ", join(supported_integrators, " "))
    end

    if driver=="LLG"
        tol = 1e-6
        field = CUDA.zeros(Float, 3*n_total)
        if integrator == "Heun"
            dopri5 = ModifiedEulerGPU(n_total, llg_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(n_total, llg_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(n_total, llg_call_back_gpu, tol)
        else
            dopri5 = DormandPrinceCayleyGPU(n_total, llg_cayley_call_back_gpu, tol)
        end
        return LLG_GPU(true, Float(0.1), Float(2.21e5), dopri5, tol, field)
    elseif driver=="LLG_STT"
        tol = 1e-6
        field = CUDA.zeros(Float, 3*n_total)
        ux = CUDA.zeros(Float, n_total)
        uy = CUDA.zeros(Float, n_total)
        uz = CUDA.zeros(Float, n_total)
        hstt = CUDA.zeros(Float, 3*n_total)
        fun = t::Float64 -> 1.0
        if integrator == "Heun"
            dopri5 = ModifiedEulerGPU(n_total, llg_stt_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(n_total, llg_stt_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(n_total, llg_stt_call_back_gpu, tol)
        else
            dopri5 = DormandPrinceCayleyGPU(n_total, llg_stt_cayley_call_back_gpu, tol)
        end
        return LLG_STT_GPU(Float(0.5), Float(0), Float(2.21e5), dopri5, tol, ux, uy, uz, hstt, field, fun)
    elseif driver=="LLG_CPP"
        tol = 1e-6
        field = CUDA.zeros(Float, 3*n_total)
        aj = CUDA.zeros(Float, n_total)
        fun = t::Float64 -> 1.0
        p = (0,0,1)
        if integrator == "Heun"
            dopri5 = ModifiedEulerGPU(n_total, llg_cpp_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(n_total, llg_cpp_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(n_total, llg_cpp_call_back_gpu, tol)
        else
            error(@sprintf("Unsupported combination driver %s and %s.", driver, integrator))
        end
        return LLG_CPP_GPU(Float(0.5), Float(0), Float(2.21e5), aj, Float(0), dopri5, tol, p, field)
    elseif driver=="LLG_STT_CPP"
        tol = 1e-6
        ux = CUDA.zeros(Float, n_total)
        uy = CUDA.zeros(Float, n_total)
        uz = CUDA.zeros(Float, n_total)
        aj = CUDA.zeros(Float, n_total)
        hstt = CUDA.zeros(Float, 3*n_total)
        field = CUDA.zeros(Float, 3*n_total)
        fun = t::Float64 -> 1.0
        p = (0,0,1)
        if integrator == "Heun"
            dopri5 = ModifiedEuler(n_total, llg_stt_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(n_total, llg_stt_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(n_total, llg_stt_cpp_call_back_gpu, tol)
        else
            dopri5 = DormandPrinceCayleyGPU(n_total, llg_stt_cpp_cayley_call_back_gpu, tol)
        end
        return LLG_STT_CPP_GPU(Float(0.5), Float(0), Float(2.21e5), Float(0), dopri5, tol, p, aj, ux, uy, uz, hstt, field, fun)
    end
    return nothing
end
