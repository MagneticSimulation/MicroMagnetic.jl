module SpinDynamics
using Printf

export create_mesh, create_sim, init_m0, add_zeeman, add_exch, add_anis, run_until, relax

include("head.jl")

include("mesh.jl")
include("llg.jl")

include("helper.jl")
include("ode.jl")

function create_sim(mesh::Mesh; name="dyn", tol=1e-6)
  nxyz = mesh.nx*mesh.ny*mesh.nz
  spin = zeros(Float64,3*nxyz)
  prespin = zeros(Float64,3*nxyz)
  field = zeros(Float64,3*nxyz)
  energy = zeros(Float64,nxyz)
  Ms = zeros(Float64,nxyz)
  dopri5 = init_runge_kutta(nxyz, rhs_call_back, tol)
  interactions = []
  return SimData(mesh, dopri5, spin, prespin, field, energy, Ms, nxyz, name, 0.1, 2.21e5, true, interactions)
end

function init_m0(sim::SimData, m0::Any; norm=true)
  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
    normalise(sim.prespin, sim.nxyz)
  end
  sim.spin[:] .= sim.prespin[:]
end

function add_zeeman(sim::SimData, H0::Any; name="zeeman")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  init_vector!(field, sim.mesh, H0)
  b = reshape(field, 3, nxyz)
  Hx = sum(b[1,:])/nxyz
  Hy = sum(b[2,:])/nxyz
  Hz = sum(b[3,:])/nxyz
  zeeman =  Zeeman(Hx, Hy, Hz, field, energy, name)
  push!(sim.interactions, zeeman)
end

function add_exch(sim::SimData, A::Float64; name="exch")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  exch =  Exchange(A, field, energy, name)
  push!(sim.interactions, exch)
end

function add_anis(sim::SimData, Ku::Float64; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Kus =  zeros(Float64, nxyz)
  Kus[:] .= Ku
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  anis =  Anisotropy(Kus, axis, field, energy, name)
  push!(sim.interactions, anis)
end

function relax(sim::SimData; maxsteps=10000, init_step = 1e-13, stopping_dmdt=0.01)
  step = 0
  sim.precession = false
  rk_data = sim.ode
  rk_data.step_next = compute_init_step(sim, init_step)
  dmdt_factor = (2 * pi / 360) * 1e9
  for i=1:maxsteps
    advance_step(sim, rk_data)
    step_size = rk_data.step
    omega_to_spin(rk_data.omega, sim.prespin, sim.spin, sim.nxyz)
    max_dmdt = compute_dmdt(sim.prespin, sim.spin, sim.nxyz, step_size)
    output = @sprintf("step =%5d   step_size=%6g    sim.t=%6g    max_dmdt=%6g", i, step_size, rk_data.t, max_dmdt/dmdt_factor)
    println(output)
    if max_dmdt < stopping_dmdt*dmdt_factor
      break
    end
  end

end

function rhs_call_back(sim::SimData, t::Float64, omega::Array{Float64})

  dw_dt = sim.ode.dw_dt
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)
  llg_rhs(dw_dt, sim.spin, sim.field, omega, sim.alpha, sim.gamma, sim.precession, sim.nxyz)

	return dw_dt

end

end
