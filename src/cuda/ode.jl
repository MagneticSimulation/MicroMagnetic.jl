function init_runge_kutta_gpu(nxyz::Int64, rhs_fun, tol::Float64)
  Float = _cuda_using_double.x ? Float64 : Float32
  omega = cuzeros(Float,3*nxyz)
  omega_t = cuzeros(Float,3*nxyz)
  dw_dt = cuzeros(Float,3*nxyz)
  k1 = cuzeros(Float,3*nxyz)
  k2 = cuzeros(Float,3*nxyz)
  k3 = cuzeros(Float,3*nxyz)
  k4 = cuzeros(Float,3*nxyz)
  k5 = cuzeros(Float,3*nxyz)
  k6 = cuzeros(Float,3*nxyz)
  k7 = cuzeros(Float,3*nxyz)
  facmax = 5.0
  facmin = 0.2
  safety = 0.824
  return Dopri5GPU(tol, 0.0, 0.0, 0.0, facmax, facmin, safety, 0, 0, omega, omega_t, dw_dt,
                   k1, k2, k3, k4, k5, k6, k7, rhs_fun, false)
end
