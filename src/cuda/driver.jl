function create_driver_gpu(driver::String, nxyz::Int64)
    if driver=="LLG"
		tol = 1e-6
        dopri5 = init_runge_kutta_gpu(nxyz, llg_call_back_gpu, tol)
        field = cuzeros(FloatGPU, 3*nxyz)
		return LLG_GPU(true, 0.1, 2.21e5, dopri5, tol, field)
    else
       error("Supported drivers for GPU: LLG.")
    end
    return nothing
end
