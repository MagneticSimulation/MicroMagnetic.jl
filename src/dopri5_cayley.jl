
mutable struct DormandPrinceCayley <: IntegratorCayley
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::Array{Float64, 1}
   omega_t::Array{Float64, 1}
   dw_dt::Array{Float64, 1}
   k1::Array{Float64, 1}
   k2::Array{Float64, 1}
   k3::Array{Float64, 1}
   k4::Array{Float64, 1}
   k5::Array{Float64, 1}
   k6::Array{Float64, 1}
   k7::Array{Float64, 1}
   rhs_fun::Function
   succeed::Bool
end


function DormandPrinceCayley(nxyz::Int64, rhs_fun, tol::Float64)
  omega = zeros(Float64,3*nxyz)
  omega_t = zeros(Float64,3*nxyz)
  dw_dt = zeros(Float64,3*nxyz)
  k1 = zeros(Float64, 3*nxyz)
  k2 = zeros(Float64, 3*nxyz)
  k3 = zeros(Float64, 3*nxyz)
  k4 = zeros(Float64, 3*nxyz)
  k5 = zeros(Float64, 3*nxyz)
  k6 = zeros(Float64, 3*nxyz)
  k7 = zeros(Float64, 3*nxyz)
  facmax = 5.0
  facmin = 0.2
  safety = 0.824
  return DormandPrinceCayley(tol, 0.0, 0, 0, facmax, facmin, safety, 0, 0, omega, omega_t, dw_dt,
                k1, k2, k3, k4, k5, k6, k7, rhs_fun, false)
end

#https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Dormand%E2%80%93Prince
function dopri5_step(sim::AbstractSim, step::Float64, t::Float64)

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

function compute_init_step(sim::AbstractSim, dt::Float64)
  abs_step = dt
  abs_step_tmp = dt
  integrator = sim.driver.ode
  fill!(integrator.omega, 0)
  integrator.rhs_fun(sim, integrator.dw_dt, integrator.t, integrator.omega)
  r_step = maximum(abs.(integrator.dw_dt)/(integrator.safety*integrator.tol^0.2))
  integrator.nfevals += 1
  #FIXME: how to obtain a reasonable init step?
  if abs_step*r_step > 0.001
    abs_step_tmp = 0.001/r_step
  end
  return min(abs_step, abs_step_tmp)
end

function advance_step(sim::AbstractSim, integrator::IntegratorCayley)

    if integrator.succeed
        omega_to_spin(integrator.omega, sim.prespin, sim.spin, sim.nxyz)
        if integrator.nsteps%10 == 0
          normalise(sim.spin, sim.nxyz)
        end
        sim.prespin .= sim.spin
    end

    t = integrator.t

    if integrator.step_next <= 0
        integrator.step_next = compute_init_step(sim, 1.0)
    end

    step_next = integrator.step_next

    while true
        max_error = dopri5_step(sim, step_next, t)/integrator.tol

        integrator.succeed = (max_error <= 1)

        if integrator.succeed
            integrator.nsteps += 1
            integrator.step = step_next
            integrator.t += integrator.step
            factor =  integrator.safety*(1.0/max_error)^0.2
            integrator.step_next = step_next*min(integrator.facmax, max(integrator.facmin, factor))
            break
        else
            factor =  integrator.safety*(1.0/max_error)^0.25
            step_next = step_next*min(integrator.facmax, max(integrator.facmin, factor))
        end
    end
  omega_to_spin(integrator.omega, sim.prespin, sim.spin, sim.nxyz)
  return true
end

function interpolation_dopri5(integrator::DormandPrinceCayley, t::Float64)
    x = (t-integrator.t+integrator.step)/integrator.step
    #assert x>=0 && x<=1
    k1,k2,k3,k4,k5,k6,k7 = integrator.k1, integrator.k2, integrator.k3, integrator.k4, integrator.k5, integrator.k6, integrator.k7
    if x == 1.0
        integrator.omega_t .= integrator.omega
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

    integrator.omega_t .= (b1.*k1 .+ b3.*k3 .+ b4.*k4 .+ b5.*k5 .+ b6.*k6 .+ b7*k7).*integrator.step
end
