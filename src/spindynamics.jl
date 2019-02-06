module SpinDynamics

include("mesh.jl")
include("llg.jl")
include("helper.jl")


mutable struct Dopri5
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::Array{Float64}
   omega_t::Array{Float64}
   dw_dt::Array{Float64}
   ks::Array{Float64}
   rhs_fun
   succeed
end


mutable struct SimData
  ode::Dopri5
  spin::Array{Float64}
  prespin::Array{Float64}
  field::Array{Float64}
  energy::Array{Float64}
  nxyz::Int64
  name::String
  alpha::Float64
  gamma::Float64
  J::Float64
  K::Float64
  Hx::Float64
  Hy::Float64
  Hz::Float64
end


function create_sim(mesh::Mesh; name="dyn", tol=1e-6)
  nxyz = mesh.nx*mesh.ny*mesh.nz
  spin = zeros(Float64,3*nxyz)
  prespin = zeros(Float64,3*nxyz)
  field = zeros(Float64,3*nxyz)
  energy = zeros(Float64,nxyz)
  dopri5 = init_runge_kutta(nxyz, rhs_call_back, tol)
  return SimData(dopri5, spin, prespin, field, energy, nxyz, name, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function init_m0(sim::SimData, mx::Float64, my::Float64, mz::Float64)
  b = reshape(sim.prespin, 3, sim.nxyz)
	for i=1:sim.nxyz
	 b[1, i] = mx
   b[2, i] = my
   b[3, i] = mz
	end
  sim.spin[:] = sim.prespin
end

function rhs_call_back(sim::SimData, t::Float64, omega::Array{Float64})

  dw_dt = sim.ode.dw_dt
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)
  llg_rhs(dw_dt, sim.spin, sim.field, omega, sim.alpha, sim.gamma, sim.nxyz)

	return dw_dt

end


include("ode.jl")

end
