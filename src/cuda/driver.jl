mutable struct EnergyMinimization_GPU{T<:AbstractFloat} <: DriverGPU
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

function create_driver_gpu(driver::String, integrator::String, nxyz::Int64)
    supported_drivers = ["SD", "LLG", "LLG_STT", "LLG_CPP", "LLG_STT_CPP"]
    if !(driver in supported_drivers)
        error("Supported drivers for GPU: ", join(supported_drivers, " "))
    end

    Float = _cuda_using_double.x ? Float64 : Float32

    if driver == "SD"
        gk = CuArrays.zeros(Float,3*nxyz)
        field = CuArrays.zeros(Float,3*nxyz)
        ss = CuArrays.zeros(Float, nxyz)
        sf = CuArrays.zeros(Float, nxyz)
        ff = CuArrays.zeros(Float, nxyz)
        return EnergyMinimization_GPU(gk, ss, sf, ff, field, Float(0.0), Float(1e-2), Float(1e-10), 0)
    end

    supported_integrators = ["Heun", "RungeKutta", "DormandPrince", "DormandPrinceCayley"]
    if !(integrator in supported_integrators)
        error("Supported integrators for GPU: ", join(supported_integrators, " "))
    end

    if driver=="LLG"
        tol = 1e-6
        field = CuArrays.zeros(Float, 3*nxyz)
        if integrator == "Heun"
            dopri5 = ModifiedEulerGPU(nxyz, llg_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(nxyz, llg_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(nxyz, llg_call_back_gpu, tol)
        else
            dopri5 = DormandPrinceCayleyGPU(nxyz, llg_cayley_call_back_gpu, tol)
        end
        return LLG_GPU(true, Float(0.1), Float(2.21e5), dopri5, tol, field)
    elseif driver=="LLG_STT"
        tol = 1e-6
        field = CuArrays.zeros(Float, 3*nxyz)
        ux = CuArrays.zeros(Float, nxyz)
        uy = CuArrays.zeros(Float, nxyz)
        uz = CuArrays.zeros(Float, nxyz)
        hstt = CuArrays.zeros(Float, 3*nxyz)
        fun = t::Float64 -> 1.0
        if integrator == "Heun"
            dopri5 = ModifiedEulerGPU(nxyz, llg_stt_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(nxyz, llg_stt_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(nxyz, llg_stt_call_back_gpu, tol)
        else
            dopri5 = DormandPrinceCayleyGPU(nxyz, llg_stt_cayley_call_back_gpu, tol)
        end
        return LLG_STT_GPU(Float(0.5), Float(0), Float(2.21e5), dopri5, tol, ux, uy, uz, hstt, field, fun)
    elseif driver=="LLG_CPP"
        tol = 1e-6
        field = CuArrays.zeros(Float, 3*nxyz)
        aj = CuArrays.zeros(Float, nxyz)
        fun = t::Float64 -> 1.0
        p = (0,0,1)
        if integrator == "Heun"
            dopri5 = ModifiedEulerGPU(nxyz, llg_cpp_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(nxyz, llg_cpp_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(nxyz, llg_cpp_call_back_gpu, tol)
        else
            error(@sprintf("Unsupported combination driver %s and %s.", driver, integrator))
        end
        return LLG_CPP_GPU(Float(0.5), Float(0), Float(2.21e5), aj, Float(0), dopri5, tol, p, field)
    elseif driver=="LLG_STT_CPP"
        tol = 1e-6
        ux = CuArrays.zeros(Float, nxyz)
        uy = CuArrays.zeros(Float, nxyz)
        uz = CuArrays.zeros(Float, nxyz)
        aj = CuArrays.zeros(Float, nxyz)
        hstt = CuArrays.zeros(Float, 3*nxyz)
        field = CuArrays.zeros(Float, 3*nxyz)
        fun = t::Float64 -> 1.0
        p = (0,0,1)
        if integrator == "Heun"
            dopri5 = ModifiedEuler(nxyz, llg_stt_call_back_gpu, 1e-14)
        elseif integrator == "RungeKutta"
            dopri5 = RungeKuttaGPU(nxyz, llg_stt_call_back_gpu, 5e-14)
        elseif integrator == "DormandPrince"
            dopri5 = DormandPrinceGPU(nxyz, llg_stt_cpp_call_back_gpu, tol)
        else
            dopri5 = DormandPrinceCayleyGPU(nxyz, llg_stt_cpp_cayley_call_back_gpu, tol)
        end
        return LLG_STT_CPP_GPU(Float(0.5), Float(0), Float(2.21e5), Float(0), dopri5, tol, p, aj, ux, uy, uz, hstt, field, fun)
    end
    return nothing
end
