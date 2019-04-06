function Sim(mesh::MeshGPU; driver="LLG", name="dyn")
  nxyz = mesh.nx*mesh.ny*mesh.nz
  spin = cuzeros(FloatGPU, 3*nxyz)
  prespin = cuzeros(FloatGPU,3*nxyz)
  field = cuzeros(FloatGPU,3*nxyz)
  energy = cuzeros(FloatGPU,nxyz)
  Ms = cuzeros(FloatGPU, nxyz)
  driver = create_driver_gpu(driver, nxyz)

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> o.total_energy, average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
  interactions = []

  blocks, threads = CuArrays.cudims(nxyz)
  return MicroSimGPU(mesh, driver, saver, spin, prespin, field, energy, Ms, nxyz, blocks, threads, name, interactions)

end

function set_Ms(sim::MicroSimGPU, init::Any)
    Ms = zeros(FloatGPU, sim.nxyz)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.Ms, Ms)
    return true
end


function init_m0(sim::MicroSimGPU, m0::Any; norm=true)
  spin = zeros(FloatGPU, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end
  copyto!(sim.spin, spin)
  copyto!(sim.prespin, sim.spin)
  return true
end

function add_zeeman(sim::MicroSimGPU, H0::Any; name="zeeman")
  nxyz = sim.nxyz
  field = zeros(FloatGPU, 3*nxyz)
  energy = zeros(FloatGPU, nxyz)
  init_vector!(field, sim.mesh, H0)
  b = reshape(field, 3, nxyz)
  Hx = sum(b[1,:])/nxyz
  Hy = sum(b[2,:])/nxyz
  Hz = sum(b[3,:])/nxyz
  field_gpu = CuArray(field)
  zeeman =  ZeemanGPU(Hx, Hy, Hz, field, field_gpu, energy, name)
  push!(sim.interactions, zeeman)

  push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
  push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
  id = length(sim.interactions)
  fun = o::AbstractSim ->  (o.interactions[id].Hx, o.interactions[id].Hy, o.interactions[id].Hz)
  push!(sim.saver.results, fun)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
end


function add_exch(sim::MicroSimGPU, A::Float64; name="exch")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  exch = ExchangeGPU(A, field, energy, FloatGPU(0.0), name)

  push!(sim.interactions, exch)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
end

function add_anis(sim::MicroSimGPU, Ku::Float64; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Kus =  cuzeros(FloatGPU, nxyz)
  Kus[:] .= Ku
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  anis =  AnisotropyGPU(Kus, axis, field, energy, name)
  push!(sim.interactions, anis)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
end


function relax(sim::MicroSimGPU; maxsteps=10000,
	     stopping_dmdt=0.01, save_m_every = 10, save_vtk_every = -1)  #TODO: merge with CPU function?
  step = 0
  rk_data = sim.driver.ode
  rk_data.step_next = compute_init_step(sim, 1e-13)
  dmdt_factor = 1.0
  if isa(sim, MicroSimGPU)
    dmdt_factor = (2 * pi / 360) * 1e9
  end
  for i=1:maxsteps
    advance_step(sim, rk_data)
    step_size = rk_data.step
    omega_to_spin(rk_data.omega, sim.prespin, sim.spin, sim.nxyz)
    println("prespin: ",sim.prespin)
    println("spin: ",sim.spin)
    max_dmdt = compute_dmdt(sim.prespin, sim.spin, sim.nxyz, step_size)
    #max_length = error_length_m(sim.spin, sim.nxyz)
    @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
	               i, step_size, rk_data.t, max_dmdt/dmdt_factor)
    if i%save_m_every == 0
      #write_data(sim)
    end
	if save_vtk_every > 0
		if i%save_vtk_every == 0
		  save_vtk(sim, @sprintf("output_%d", i))
		end
	end
    sim.saver.t = rk_data.t
    sim.saver.nsteps += 1
    if max_dmdt < stopping_dmdt*dmdt_factor
      @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
      break
    end
  end
end
