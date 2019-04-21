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
end

function create_driver_gpu(driver::String, nxyz::Int64)
    if driver == "SD"
        Float = _cuda_using_double.x ? Float64 : Float32
		gk = cuzeros(Float,3*nxyz)
		field = cuzeros(Float,3*nxyz)
        ss = cuzeros(Float, nxyz)
        sf = cuzeros(Float, nxyz)
        ff = cuzeros(Float, nxyz)
		return EnergyMinimization_GPU(gk, ss, sf, ff, field, Float(0.0), Float(1e-2), Float(1e-10), 0)
    elseif driver=="LLG"
		tol = 1e-6
        dopri5 = init_runge_kutta_gpu(nxyz, llg_call_back_gpu, tol)
        Float = _cuda_using_double.x ? Float64 : Float32
        field = cuzeros(Float, 3*nxyz)
		return LLG_GPU(true, Float(0.1), Float(2.21e5), dopri5, tol, field)
	elseif driver=="LLG_STT"
        tol = 1e-6
		dopri5 = init_runge_kutta_gpu(nxyz, llg_stt_call_back_gpu, tol)
		Float = _cuda_using_double.x ? Float64 : Float32
        ux = cuzeros(Float, nxyz)
        uy = cuzeros(Float, nxyz)
        uz = cuzeros(Float, nxyz)
        hstt = cuzeros(Float, 3*nxyz)
        field = cuzeros(Float, 3*nxyz)
        return LLG_STT_GPU(Float(0.5), Float(0), Float(2.21e5), dopri5, tol, ux, uy, uz, hstt, field)
    else
       error("Supported drivers for GPU: LLG, LLG_STT")
    end
    return nothing
end
