include("field.jl")
#include("llg.jl")

function init_runge_kutta(nxyz::Int64, rhs_fun, tol::Float64)
  omega = zeros(Float64,3*nxyz)
  omega_t = zeros(Float64,3*nxyz)
  dw_dt = zeros(Float64,3*nxyz)
  ks = zeros(Float64, 3*nxyz, 7)
  facmax = 5.0
  facmin = 0.2
  safety = 0.824
  return Dopri5(tol, 0.0, 0, 0, facmax, facmin, safety, 0, 0, omega, omega_t, dw_dt, ks, rhs_fun, false)
end

function dopri5_step(sim::AbstractSim, step::Float64, t::Float64)
  a = (1/5, 3/10, 4/5, 8/9, 1, 1)
  b = (1/5, 3/40, 9/40, 44/45, -56/15, 32/9)
  c = (19372/6561, -25360/2187, 64448/6561, -212/729)
  d = (9017/3168, -355/33, 46732/5247, 49/176, -5103/18656)
  v = (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84)
  w = (71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40)
  ode = sim.driver.ode
  rhs = ode.rhs_fun
  ks = ode.ks
  y_next = ode.omega

  fill!(y_next, 0) # we always have y=0
  ks[:,1] = rhs(sim, t, y_next)  # we can not simply copy k7 since we have modified y each time

  y_next[:] = b[1]*ks[:,1]*step
  ks[:,2] = rhs(sim, t + a[1]*step, y_next) #k2

  y_next[:] = (b[2]*ks[:,1]+b[3]*ks[:,2])*step
  ks[:,3] = rhs(sim, t + a[2]*step, y_next) #k3

  y_next[:] = (b[4]*ks[:,1]+b[5]*ks[:,2]+b[6]*ks[:,3])*step
  ks[:,4] = rhs(sim, t + a[3]*step, y_next) #k4

  y_next[:] = (c[1]*ks[:,1]+c[2]*ks[:,2]+c[3]*ks[:,3]+c[4]*ks[:,4])*step
  ks[:,5] = rhs(sim, t + a[4]*step, y_next) #k5

  y_next[:] = (d[1]*ks[:,1]+d[2]*ks[:,2]+d[3]*ks[:,3]+d[4]*ks[:,4]+d[5]*ks[:,5])*step
  ks[:,6] = rhs(sim, t + a[5]*step, y_next) #k6

  y_next[:] = (v[1]*ks[:,1]+v[2]*ks[:,2]+v[3]*ks[:,3]+v[4]*ks[:,4]+v[5]*ks[:,5]+v[6]*ks[:,6])*step
  ks[:,7] = rhs(sim, t + a[6]*step, y_next) #k7

  ode.nfevals += 7
  error = ode.omega_t #we make use of omega_t to store the error temporary
  error[:] = (w[1]*ks[:,1]+w[2]*ks[:,2]+w[3]*ks[:,3]+w[4]*ks[:,4]+w[5]*ks[:,5]+w[6]*ks[:,6]+w[7]*ks[:,7])*step

  max_error =  maximum(abs.(error)) + eps()

  #omega_to_spin(error, sim.prespin, sim.spin, sim.nxyz)
  #sim.spin[:] -= sim.prespin[:]
  #max_dm_error = maximum(abs.(sim.spin)) + eps()

  return max_error
end

function compute_init_step(sim::AbstractSim, dt::Float64)
  abs_step = dt
  abs_step_tmp = dt
  rk_data = sim.driver.ode
  fill!(rk_data.omega, 0)
  r_step = maximum(abs.(rk_data.rhs_fun(sim, rk_data.t, rk_data.omega))/(rk_data.safety*rk_data.tol^0.2))
  rk_data.nfevals += 1
  if abs_step*r_step > 1
    abs_step_tmp = 1.0/r_step
  end
  return min(abs_step, abs_step_tmp)
end

function advance_step(sim::AbstractSim, rk_data::Dopri5)

    if rk_data.succeed
        omega_to_spin(rk_data.omega, sim.prespin, sim.spin, sim.nxyz)
        if rk_data.nsteps%10 == 0
          normalise(sim.spin, sim.nxyz)
        end
        sim.prespin[:] = sim.spin[:]
    end

    t = rk_data.t

    if rk_data.step_next <= 0
        rk_data.step_next = compute_init_step(sim, 1.0)
    end

    step_next = rk_data.step_next

    while true
        max_error = dopri5_step(sim, step_next, t)/rk_data.tol

        rk_data.succeed = (max_error <= 1)

        if rk_data.succeed
            rk_data.nsteps += 1
            rk_data.step = step_next
            rk_data.t += rk_data.step
            factor =  rk_data.safety*(1.0/max_error)^0.2
            rk_data.step_next = step_next*min(rk_data.facmax, max(rk_data.facmin, factor))
            break
        else
            factor =  rk_data.safety*(1.0/max_error)^0.25
            step_next = step_next*min(rk_data.facmax, max(rk_data.facmin, factor))
        end
    end
	omega_to_spin(rk_data.omega, sim.prespin, sim.spin, sim.nxyz)
end

function interpolation_dopri5(rk_data::Dopri5, t::Float64)
    x = (t-rk_data.t+rk_data.step)/rk_data.step
    #assert x>=0 && x<=1
    ks = rk_data.ks
    if x == 1.0
        rk_data.omega_t[:] = rk_data.omega[:]
        return
    end

    v = (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84)
    x1 = x*x*(3-2*x)
    x2 = x*x*(x-1)^2
    b1 = x1*v[1] + x*(x-1)^2 - x2*5*(2558722523 - 31403016*x)/11282082432
    b2 = 0
    b3 = x1*v[3] + x2*100*(882725551 - 15701508*x)/32700410799
    b4 = x1*v[4] - x2*25*(443332067 - 31403016*x)/1880347072
    b5 = x1*v[5] + x2*32805*(23143187 - 3489224*x)/199316789632
    b6 = x1*v[6] - x2*55*(29972135 - 7076736*x)/822651844
    b7 = x*x*(x-1) + x2*10*(7414447 - 829305*x)/29380423

    rk_data.omega_t[:] = (b1*ks[:,1]+b2*ks[:,2]+b3*ks[:,3]+b4*ks[:,4]+b5*ks[:,5]+b6*ks[:,6]+b7*ks[:,7])*rk_data.step
end
