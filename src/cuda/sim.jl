function FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc="open")
  nxyz = nx*ny*nz
  volume = dx*dy*dz
  Float = _cuda_using_double.x ? Float64 : Float32
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false
  return FDMeshGPU(Float(dx), Float(dy), Float(dz), nx, ny, nz, nxyz, xperiodic, yperiodic, zperiodic, Float(volume))
end

function Sim(mesh::MeshGPU; driver="LLG", name="dyn")
  nxyz = mesh.nx*mesh.ny*mesh.nz
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = cuzeros(Float, 3*nxyz)
  prespin = cuzeros(Float,3*nxyz)
  field = cuzeros(Float,3*nxyz)
  energy = cuzeros(Float,nxyz)
  Ms = cuzeros(Float, nxyz)
  driver = create_driver_gpu(driver, nxyz)

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> sum(o.energy), average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
  interactions = []

  blocks, threads = CuArrays.cudims(nxyz)
  return MicroSimGPU(mesh, driver, saver, spin, prespin, field, energy, Ms, Float(0.0),
                     nxyz, blocks, threads, name, interactions)

end

function set_Ms(sim::MicroSimGPU, init::Any)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms = zeros(Float, sim.nxyz)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.Ms, Ms)
    return true
end


function init_m0(sim::MicroSimGPU, m0::Any; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
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
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  init_vector!(field, sim.mesh, H0)
  b = reshape(field, 3, nxyz)
  Hx = sum(b[1,:])/nxyz
  Hy = sum(b[2,:])/nxyz
  Hz = sum(b[3,:])/nxyz
  field_gpu = CuArray(field)
  zeeman =  ZeemanGPU(Float64(Hx), Float64(Hy), Float64(Hz), field, field_gpu, energy, Float(0.0), name)
  push!(sim.interactions, zeeman)

  push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
  push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
  id = length(sim.interactions)
  fun = o::AbstractSim ->  (o.interactions[id].Hx, o.interactions[id].Hy, o.interactions[id].Hz)
  push!(sim.saver.results, fun)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return zeeman
end


function add_exch(sim::MicroSimGPU, A::Float64; name="exch")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  exch = ExchangeGPU(Float(A), field, energy, Float(0.0), name)

  push!(sim.interactions, exch)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return exch
end

function add_dmi(sim::MicroSimGPU, D::Float64; name="dmi")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  dmi = BulkDMIGPU(Float(D), field, energy, Float(0.0), name)

  push!(sim.interactions, dmi)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return dmi
end

function add_anis(sim::MicroSimGPU, Ku::Float64; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  Kus =  cuzeros(Float, nxyz)
  Kus[:] .= Ku
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  anis =  AnisotropyGPU(Kus, axis, field, energy, Float(0.0), name)
  push!(sim.interactions, anis)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return anis
end

function add_demag(sim::MicroSimGPU; name="demag")
  demag = init_demag_gpu(sim)
  demag.name = name
  push!(sim.interactions, demag)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::MicroSimGPU->sum(o.interactions[id].energy))
  return demag
end

function relax(sim::MicroSimGPU; maxsteps=10000, init_step = 1e-13, stopping_dmdt=0.01, stopping_torque=0.1, save_m_every = 10, save_vtk_every=-1)
  if isa(sim.driver, EnergyMinimization_GPU)
    relax(sim, sim.driver, maxsteps=maxsteps, stopping_torque=stopping_torque, save_m_every=save_m_every, save_vtk_every=save_vtk_every)
elseif isa(sim.driver, LLG_GPU)
    relax(sim, sim.driver, maxsteps=maxsteps, stopping_dmdt=stopping_dmdt, save_m_every=save_m_every, save_vtk_every=save_vtk_every)
  end
  return nothing
end

function relax(sim::MicroSimGPU, driver::EnergyMinimization_GPU; maxsteps=10000,
               stopping_torque=0.1, save_m_every = 10, save_vtk_every = -1)
  T = _cuda_using_double.x ? Float64 : Float32
  gk_abs = cuzeros(T, 3*sim.nxyz)
  for i=1:maxsteps
    run_step(sim, sim.driver)
    abs!(gk_abs, driver.gk)  #max_torque = maximum(abs.(driver.gk)) eats gpu memory???
    max_torque = maximum(gk_abs)
	@info @sprintf("step=%5d  tau=%10.6e  max_torque=%10.6e", driver.steps, driver.tau, max_torque)
    if i%save_m_every == 0
      compute_system_energy(sim, sim.spin, 0.0)
      write_data(sim)
    end
	if save_vtk_every > 0
		if i%save_vtk_every == 0
	  	  save_vtk(sim, @sprintf("%s_%d", sim.name, i))
		end
	end
    sim.saver.nsteps += 1
    if max_torque < stopping_torque
      @info @sprintf("max_torque (mxmxH) is less than stopping_torque=%g, Done!", stopping_torque)
      break
    end
  end
end

function relax(sim::MicroSimGPU, driver::LLG_GPU; maxsteps=10000, stopping_dmdt=0.01, save_m_every = 10, save_vtk_every = -1)  #TODO: merge with CPU function?
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
    compute_dm(rk_data.omega_t, sim.prespin, sim.spin, sim.nxyz)
	max_dmdt = maximum(rk_data.omega_t)/step_size
    #max_length = error_length_m(sim.spin, sim.nxyz)
    @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
	               i, step_size, rk_data.t, max_dmdt/dmdt_factor)
    if i%save_m_every == 0
      compute_system_energy(sim, sim.spin, 0.0)
      write_data(sim)
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
