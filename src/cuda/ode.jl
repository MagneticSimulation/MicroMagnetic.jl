function init_runge_kutta_gpu(nxyz::Int64, rhs_fun, tol::Float64)
  omega = cuzeros(FloatGPU,3*nxyz)
  omega_t = cuzeros(FloatGPU,3*nxyz)
  dw_dt = cuzeros(FloatGPU,3*nxyz)
  k1 = cuzeros(FloatGPU,3*nxyz)
  k2 = cuzeros(FloatGPU,3*nxyz)
  k3 = cuzeros(FloatGPU,3*nxyz)
  k4 = cuzeros(FloatGPU,3*nxyz)
  k5 = cuzeros(FloatGPU,3*nxyz)
  k6 = cuzeros(FloatGPU,3*nxyz)
  k7 = cuzeros(FloatGPU,3*nxyz)
  facmax = 5.0
  facmin = 0.2
  safety = 0.824
  return Dopri5GPU(tol, 0.0, 0, 0, facmax, facmin, safety, 0, 0, omega, omega_t, dw_dt, k1, k2, k3, k4, k5, k6, k7, rhs_fun, false)
end

#https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Dormand%E2%80%93Prince
function dopri5_step(sim::MicroSimGPU, step::Float64, t::Float64)

  a = (1/5, 3/10, 4/5, 8/9, 1.0, 1.0)
  b = (1/5, 3/40, 9/40, 44/45, -56/15, 32/9)
  c = (19372/6561, -25360/2187, 64448/6561, -212/729)
  d = (9017/3168, -355/33, 46732/5247, 49/176, -5103/18656)
  v = (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84)
  w = (71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40)
  ode = sim.driver.ode
  rhs = ode.rhs_fun
  y_next = ode.omega
  k1,k2,k3,k4,k5,k6,k7 = ode.k1, ode.k2, ode.k3, ode.k4, ode.k5, ode.k6, ode.k7

  fill!(y_next, 0) # we always have y=0
  ode.rhs_fun(sim, k1, t, y_next) #compute k1

  y_next .= b[1].*k1.*step
  ode.rhs_fun(sim, k2, t + a[1]*step, y_next) #k2

  y_next .= (b[2].*k1 .+ b[3].*k2).*step
  ode.rhs_fun(sim, k3, t + a[2]*step, y_next) #k3

  y_next .= (b[4].*k1 .+ b[5].*k2 .+ b[6].*k3).*step
  ode.rhs_fun(sim, k4, t + a[3]*step, y_next) #k4

  y_next .= (c[1].*k1 .+ c[2].*k2 + c[3].*k3 .+ c[4].*k4).*step
  ode.rhs_fun(sim, k5, t + a[4]*step, y_next) #k5

  y_next .= (d[1].*k1 .+ d[2].*k2 .+ d[3].*k3 .+ d[4].*k4 + d[5].*k5).*step
  ode.rhs_fun(sim, k6, t + a[5]*step, y_next) #k6

  y_next .= (v[1].*k1 .+ v[2].*k2 .+ v[3].*k3 .+ v[4].*k4 .+ v[5].*k5 .+ v[6].*k6) .* step
  ode.rhs_fun(sim, k7, t + a[6]*step, y_next) #k7

  ode.nfevals += 7
  error = ode.omega_t #we make use of omega_t to store the error temporary
  error .= (w[1].*k1 + w[2].*k2 .+ w[3].*k3 .+ w[4].*k4 .+ w[5].*k5 + w[6].*k6 + w[7].*k7).*step

  max_error =  maximum(abs.(error)) + eps()

  return max_error
end

function compute_init_step(sim::MicroSimGPU, dt::Float64)
  abs_step = dt
  abs_step_tmp = dt
  rk_data = sim.driver.ode
  fill!(rk_data.omega, 0)
  rk_data.rhs_fun(sim, rk_data.dw_dt, rk_data.t, rk_data.omega)
  r_step = maximum(abs.(rk_data.dw_dt)/(rk_data.safety*rk_data.tol^0.2))
  rk_data.nfevals += 1
  if abs_step*r_step > 1
    abs_step_tmp = 1.0/r_step
  end
  return min(abs_step, abs_step_tmp)
end

function advance_step(sim::MicroSimGPU, rk_data::Dopri5GPU)
    if rk_data.succeed
        omega_to_spin(rk_data.omega, sim.prespin, sim.spin, sim.nxyz)
        if rk_data.nsteps%10 == 0
          normalise(sim.spin, sim.nxyz)
        end
        sim.prespin .= sim.spin
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
end
