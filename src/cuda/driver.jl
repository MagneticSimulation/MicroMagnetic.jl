function create_driver_gpu(driver::String, nxyz::Int64)
    if driver=="LLG"
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
