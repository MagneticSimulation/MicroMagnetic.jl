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

#y = c1*k1*step
function rk__step1__kernel!(y, c1, k1, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
       @inbounds y[i] = c1*k1[i]*step
    end
	return nothing
end

#y = (c1*k1+c2*k2)*step
function rk__step2__kernel!(y, c1, k1, c2, k2, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i])*step
    end
	return nothing
end

function rk__step3__kernel!(y, c1, k1, c2, k2, c3, k3, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i]+c3*k3[i])*step
    end
	return nothing
end

function rk__step4__kernel!(y, c1, k1, c2, k2, c3, k3, c4, k4, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i])*step
    end
	return nothing
end

function rk__step5__kernel!(y, c1, k1, c2, k2, c3, k3, c4, k4, c5, k5, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i])*step
    end
	return nothing
end

function rk__step6__kernel!(y, c1, k1, c2, k2, c3, k3, c4, k4, c5, k5, c6, k6, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i])*step
    end
	return nothing
end

function rk__step7__kernel!(y, c1, k1, c2, k2, c3, k3, c4, k4, c5, k5, c6, k6, c7, k7, step, n)
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]+c7*k7[i])*step
    end
	return nothing
end

function dopri5_step(sim::MicroSimGPU, step::Float64, t::Float64)

  a = (1/5, 3/10, 4/5, 8/9, 1.0, 1.0)
  b = (1/5, 3/40, 9/40, 44/45, -56/15, 32/9)
  c = (19372/6561, -25360/2187, 64448/6561, -212/729)
  d = (9017/3168, -355/33, 46732/5247, 49/176, -5103/18656)
  v = (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84)
  w = (71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40)
  ode = sim.driver.ode

  N = length(ode.omega)
  blk, thr = CuArrays.cudims(ode.omega)

  fill!(ode.omega, 0) # we always have y=0
  ode.rhs_fun(sim, ode.k1, t, ode.omega) #compute k1

  @cuda blocks=blk threads=thr rk__step1__kernel!(ode.omega,
                               b[1], ode.k1, step, N)
  ode.rhs_fun(sim, ode.k2, t + a[1]*step, ode.omega) #k2

  @cuda blocks=blk threads=thr rk__step2__kernel!(ode.omega,
                               b[2], ode.k1, b[3], ode.k2,
                               step, N)
  ode.rhs_fun(sim, ode.k3, t + a[2]*step, ode.omega) #k3

  @cuda blocks=blk threads=thr rk__step3__kernel!(ode.omega,
                               b[4], ode.k1, b[5], ode.k2,
                               b[6], ode.k3, step, N)
  ode.rhs_fun(sim, ode.k4, t + a[3]*step, ode.omega) #k4

  @cuda blocks=blk threads=thr rk__step4__kernel!(ode.omega,
                               c[1], ode.k1, c[2], ode.k2,
                               c[3], ode.k3, c[4], ode.k4,
                               step, N)
  ode.rhs_fun(sim, ode.k5, t + a[4]*step, ode.omega) #k5

  @cuda blocks=blk threads=thr rk__step5__kernel!(ode.omega,
                               d[1], ode.k1, d[2], ode.k2,
                               d[3], ode.k3, d[4], ode.k4,
                               d[5], ode.k5, step, N)
  ode.rhs_fun(sim, ode.k6, t + a[5]*step, ode.omega) #k6

  @cuda blocks=blk threads=thr rk__step6__kernel!(ode.omega,
                               v[1], ode.k1, v[2], ode.k2,
                               v[3], ode.k3, v[4], ode.k4,
                               v[5], ode.k5, v[6], ode.k6,
                               step, N)
  ode.rhs_fun(sim, ode.k7, t + a[6]*step, ode.omega) #k7

  ode.nfevals += 7
  error = ode.omega_t #we make use of omega_t to store the error temporary
  @cuda blocks=blk threads=thr rk__step7__kernel!(error,
                               w[1], ode.k1, w[2], ode.k2,
                               w[3], ode.k3, w[4], ode.k4,
                               w[5], ode.k5, w[6], ode.k6,
                               w[7], ode.k5, step, N)
  abs!(error)
  max_error =  maximum(error) + eps()

  return max_error
end

function compute_init_step(sim::MicroSimGPU, dt::Float64)
  abs_step = dt
  abs_step_tmp = dt
  rk_data = sim.driver.ode
  fill!(rk_data.omega, 0)
  rk_data.rhs_fun(sim, rk_data.dw_dt, rk_data.t, rk_data.omega)
  abs!(rk_data.dw_dt)
  r_step = maximum(rk_data.dw_dt)/(rk_data.safety*rk_data.tol^0.2)
  rk_data.nfevals += 1
  if abs_step*r_step > 1
    abs_step_tmp = 1.0/r_step
  end
  return min(abs_step, abs_step_tmp)
end


function interpolation_dopri5(ode::Dopri5GPU, t::Float64)
    x = (t-ode.t+ode.step)/ode.step
    #assert x>=0 && x<=1
    if x == 1.0
        ode.omega_t .= ode.omega
        return
    end

    v = (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84)
    x1 = x*x*(3-2*x)
    x2 = x*x*(x-1)^2
    b1 = x1*v[1] + x*(x-1)^2 - x2*5*(2558722523 - 31403016*x)/11282082432
    #b2 = 0
    b3 = x1*v[3] + x2*100*(882725551 - 15701508*x)/32700410799
    b4 = x1*v[4] - x2*25*(443332067 - 31403016*x)/1880347072
    b5 = x1*v[5] + x2*32805*(23143187 - 3489224*x)/199316789632
    b6 = x1*v[6] - x2*55*(29972135 - 7076736*x)/822651844
    b7 = x*x*(x-1) + x2*10*(7414447 - 829305*x)/29380423

    N = length(ode.omega)
    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr rk__step6__kernel!(ode.omega_t,
                                 b1, ode.k1, b3, ode.k3,
                                 b4, ode.k4, b5, ode.k5,
                                 b6, ode.k6, b7, ode.k7,
                                 ode.step, N)
    return nothing
end
