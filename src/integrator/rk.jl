#implement the classical Runge-Kutta method for testing purpose

mutable struct RungeKutta <: Integrator
   t::Float64
   step::Float64
   nsteps::Int64
   k1::Array{Float64, 1}
   k2::Array{Float64, 1}
   k3::Array{Float64, 1}
   k4::Array{Float64, 1}
   rhs_fun::Function
end

function RungeKutta(nxyz::Int64, rhs_fun, step::Float64)
  k1 = zeros(Float64,3*nxyz)
  k2 = zeros(Float64,3*nxyz)
  k3 = zeros(Float64,3*nxyz)
  k4 = zeros(Float64,3*nxyz)
  return RungeKutta(0.0, step, 0, k1, k2, k3, k4, rhs_fun)
end

function advance_step(sim::AbstractSim, rk_data::RungeKutta)
    h = rk_data.step
	k1 =  rk_data.k1
	k2 =  rk_data.k2
	k3 =  rk_data.k3
	k4 =  rk_data.k4

	#compute k1
    rk_data.rhs_fun(sim, k1, sim.spin, rk_data.t)

    #compute k2
	sim.prespin .= sim.spin .+ 0.5 .*h.*k1
	rk_data.rhs_fun(sim, k2, sim.prespin, rk_data.t+0.5*h)

    #compute k3
	sim.prespin .= sim.spin .+ 0.5 .*h.*k2
	rk_data.rhs_fun(sim, k3, sim.prespin, rk_data.t+0.5*h)

    #compute k4
	sim.prespin .= sim.spin .+ h.*k3
	rk_data.rhs_fun(sim, k4, sim.prespin, rk_data.t+h)

    sim.prespin .= sim.spin
    sim.spin .+= (1.0/6).*(k1 .+ 2 .*k2 + 2 .*k3 .+ k4).*h

	normalise(sim.spin, sim.nxyz)
	rk_data.t += h

end
