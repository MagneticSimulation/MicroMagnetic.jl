function init_runge_kutta_gpu(nxyz::Int64, rhs_fun, tol::Float64)
  omega = cuzeros(FloatGPU,3*nxyz)
  omega_t = cuzeros(FloatGPU,3*nxyz)
  dw_dt = cuzeros(FloatGPU,3*nxyz)
  ks = cuzeros(FloatGPU, 3*nxyz, 7)
  facmax = 5.0
  facmin = 0.2
  safety = 0.824
  return Dopri5GPU(tol, 0.0, 0, 0, facmax, facmin, safety, 0, 0, omega, omega_t, dw_dt, ks, rhs_fun, false)
end
