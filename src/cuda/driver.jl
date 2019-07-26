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
  ode::Dopri5GPU
  tol::Float64
  field::CuArray{T, 1}
end

mutable struct LLG_STT_GPU{T<:AbstractFloat} <: DriverGPU
  alpha::T
  beta::T
  gamma::T
  ode::Dopri5GPU
  tol::Float64
  ux::CuArray{T, 1}
  uy::CuArray{T, 1}
  uz::CuArray{T, 1}
  h_stt::CuArray{T, 1}
  field::CuArray{T, 1}
  ufun::Function
end

mutable struct LLG_STT_CPP_GPU{T<:AbstractFloat} <: DriverGPU
  alpha::T
  beta::T
  gamma::T
  eta::T
  ode::Dopri5GPU
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

function create_driver_gpu(driver::String, nxyz::Int64)
    if driver == "SD"
        Float = _cuda_using_double.x ? Float64 : Float32
		gk = CuArrays.zeros(Float,3*nxyz)
		field = CuArrays.zeros(Float,3*nxyz)
        ss = CuArrays.zeros(Float, nxyz)
        sf = CuArrays.zeros(Float, nxyz)
        ff = CuArrays.zeros(Float, nxyz)
		return EnergyMinimization_GPU(gk, ss, sf, ff, field, Float(0.0), Float(1e-2), Float(1e-10), 0)
    elseif driver=="LLG"
		tol = 1e-6
        dopri5 = init_runge_kutta_gpu(nxyz, llg_call_back_gpu, tol)
        Float = _cuda_using_double.x ? Float64 : Float32
        field = CuArrays.zeros(Float, 3*nxyz)
		return LLG_GPU(true, Float(0.1), Float(2.21e5), dopri5, tol, field)
	elseif driver=="LLG_STT"
        tol = 1e-6
        dopri5 = init_runge_kutta_gpu(nxyz, llg_stt_call_back_gpu, tol)
        Float = _cuda_using_double.x ? Float64 : Float32
        ux = CuArrays.zeros(Float, nxyz)
        uy = CuArrays.zeros(Float, nxyz)
        uz = CuArrays.zeros(Float, nxyz)
        hstt = CuArrays.zeros(Float, 3*nxyz)
        field = CuArrays.zeros(Float, 3*nxyz)
        fun = t::Float64 -> 1.0
        return LLG_STT_GPU(Float(0.5), Float(0), Float(2.21e5), dopri5, tol, ux, uy, uz, hstt, field, fun)
    elseif driver=="LLG_STT_CPP"
        tol = 1e-6
		dopri5 = init_runge_kutta_gpu(nxyz, llg_stt_cpp_call_back_gpu, tol)
		Float = _cuda_using_double.x ? Float64 : Float32
        ux = CuArrays.zeros(Float, nxyz)
        uy = CuArrays.zeros(Float, nxyz)
        uz = CuArrays.zeros(Float, nxyz)
        aj = CuArrays.zeros(Float, nxyz)
        hstt = CuArrays.zeros(Float, 3*nxyz)
        field = CuArrays.zeros(Float, 3*nxyz)
        fun = t::Float64 -> 1.0
        p = (0,0,1)
        return LLG_STT_CPP_GPU(Float(0.5), Float(0), Float(2.21e5), Float(0), dopri5, tol, p, aj, ux, uy, uz, hstt, field, fun)
    else
       error("Supported drivers for GPU: LLG, LLG_STT, LLG_STT_CPP")
    end
    return nothing
end
