function FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc=(false, false, false))
  nxyz = nx*ny*nz
  volume = dx*dy*dz
  Float = _cuda_using_double.x ? Float64 : Float32
  return FDMeshGPU(Float(dx), Float(dy), Float(dz), nx, ny, nz, nxyz, pbc[1], pbc[2], pbc[3], Float(volume))
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
  return MicroSimGPU(mesh, driver, saver, spin, prespin, field, energy, Ms, nxyz, blocks, threads, name, interactions)

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
end

function add_demag(sim::MicroSimGPU; name="demag")
  demag = init_demag_gpu(sim)
  demag.name = name
  push!(sim.interactions, demag)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::MicroSim->sum(o.interactions[id].energy))
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
    max_dmdt = compute_dmdt(sim.prespin, sim.spin, sim.nxyz, step_size)
    #max_length = error_length_m(sim.spin, sim.nxyz)
    @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
	               i, step_size, rk_data.t, max_dmdt/dmdt_factor)
    if i%save_m_every == 0
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
