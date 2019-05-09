function Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="Dopri5")
  nxyz = mesh.nx*mesh.ny*mesh.nz
  spin = zeros(Float64,3*nxyz)
  prespin = zeros(Float64,3*nxyz)
  field = zeros(Float64,3*nxyz)
  energy = zeros(Float64,nxyz)
  Ms = zeros(Float64,nxyz)
  driver = create_driver(driver, integrator, nxyz)

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> sum(o.energy), average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
  interactions = []
  if isa(mesh, FDMesh)
    return MicroSim(mesh, driver, saver, spin, prespin, field, energy, Ms, nxyz, name, interactions)
  else
    return AtomicSim(mesh, driver, saver, spin, prespin, field, energy, Ms, nxyz, name, interactions)
  end
end

function set_Ms(sim::MicroSim, fun_Ms::Function)
    mesh = sim.mesh
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        sim.Ms[id] = fun_Ms(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    return true
end

function set_Ms(sim::MicroSim, Ms::Number)
    for i =1:sim.nxyz
        sim.Ms[i] = Ms
    end
    return true
end

function set_mu_s(sim::AtomicSim, fun_Ms::Function)
    mesh = sim.mesh
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        sim.mu_s[id] = fun_Ms(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    return true
end

function set_mu_s(sim::AtomicSim, Ms::Number)
    for i =1:sim.nxyz
        sim.mu_s[i] = Ms
    end
    return true
end

function set_ux(sim::AbstractSim, init_ux)
	init_scalar!(sim.driver.ux, sim.mesh, init_ux)
end

function average_m(sim::AbstractSim)
  b = reshape(sim.spin, 3, sim.nxyz)
  mx,my,mz = 0.0,0.0,0.0
  n = 0
  ms = isa(sim, MicroSim) ? sim.Ms : sim.mu_s
  for i = 1:sim.nxyz
    if ms[i]>0
      n += 1
      mx += b[1,i]
      my += b[2,i]
      mz += b[3,i]
    end
  end
  if n == 0
    error("n should not be zero!")
  end
  return (mx/n, my/n, mz/n)
end

function init_m0(sim::AbstractSim, m0::Any; norm=true)
  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
    normalise(sim.prespin, sim.nxyz)
  end
  sim.spin[:] .= sim.prespin[:]
end

function add_zeeman(sim::AbstractSim, H0::Any; name="zeeman")
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

  push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
  push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
  id = length(sim.interactions)
  fun = o::AbstractSim ->  (o.interactions[id].Hx, o.interactions[id].Hy, o.interactions[id].Hz)
  push!(sim.saver.results, fun)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  return zeeman
end

function add_exch(sim::AbstractSim, A::Float64; name="exch")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  if isa(sim, MicroSim)
    exch = Exchange(A, field, energy, name)
  else
	exch = HeisenbergExchange(A, field, energy, name)
  end
  push!(sim.interactions, exch)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  return exch
end

function add_dmi(sim::AbstractSim, D::Tuple{Real, Real, Real}; name="dmi")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  dmi =  BulkDMI(Float64(D[1]), Float64(D[2]), Float64(D[3]), field, energy, name)
  push!(sim.interactions, dmi)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  return dmi
end

function add_dmi(sim::AbstractSim, D::Real; name="dmi", type="bulk")
    if type == "interfacial"
        return add_dmi_interfacial(sim, D, name=name)
    end
   return add_dmi(sim, (D,D,D), name=name)
end

function add_dmi_interfacial(sim::AbstractSim, D::Real; name="dmi")
    nxyz = sim.nxyz
    field = zeros(Float64, 3*nxyz)
    energy = zeros(Float64, nxyz)
    dmi =  InterfacialDMI(Float64(D), field, energy, name)
    push!(sim.interactions, dmi)

    push!(sim.saver.headers, string("E_",name))
    push!(sim.saver.units, "J")
    id = length(sim.interactions)
    push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
    return dmi
end

function add_demag(sim::MicroSim; name="demag")
  demag = init_demag(sim)
  demag.name = name
  push!(sim.interactions, demag)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::MicroSim->sum(o.interactions[id].energy))
  return demag
end


function add_anis(sim::AbstractSim, init_ku::Any; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Kus =  zeros(Float64, nxyz)
  init_scalar!(Kus, sim.mesh, init_ku)
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  anis =  Anisotropy(Kus, axis, field, energy, name)
  push!(sim.interactions, anis)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  return anis
end

function relax(sim::AbstractSim; maxsteps=10000, init_step = 1e-13, stopping_dmdt=0.01, stopping_torque=0.1, save_m_every = 10, save_vtk_every=-1)
  if isa(sim.driver, EnergyMinimization)
    relax(sim, sim.driver, maxsteps=maxsteps, stopping_torque=stopping_torque, save_m_every=save_m_every, save_vtk_every=save_vtk_every)
  elseif isa(sim.driver, LLG)
    relax(sim, sim.driver, maxsteps=maxsteps, stopping_dmdt=stopping_dmdt, save_m_every=save_m_every, save_vtk_every=save_vtk_every)
  end
  return nothing
end

function relax(sim::AbstractSim, driver::EnergyMinimization; maxsteps=10000,
               stopping_torque=0.1, save_m_every = 10, save_vtk_every = -1)
  for i=1:maxsteps
    run_step(sim, sim.driver)
	max_torque = maximum(abs.(driver.gk))
    max_length_error = error_length_m(sim.spin, sim.nxyz)
	@info @sprintf("step=%5d  tau=%10.6e  max_torque=%10.6e  |m|_error=%10.6e",
	               i, driver.tau, max_torque, max_length_error)
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

function relax(sim::AbstractSim, driver::LLG; maxsteps=10000,
	     stopping_dmdt=0.01, save_m_every = 10, save_vtk_every = -1)
  step = 0
  rk_data = sim.driver.ode
  #rk_data.step_next = compute_init_step(sim, 1e-13)
  dmdt_factor = 1.0
  if isa(sim, MicroSim)
    dmdt_factor = (2 * pi / 360) * 1e9
  end
  for i=1:maxsteps
    advance_step(sim, rk_data)
    max_dmdt = compute_dmdt(sim.prespin, sim.spin, sim.nxyz, rk_data.step)
    max_length = error_length_m(sim.spin, sim.nxyz)
    @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e  |m|_error=%10.6e",
	               i, rk_data.step, rk_data.t, max_dmdt/dmdt_factor, max_length)
    if i%save_m_every == 0
	  compute_system_energy(sim, sim.spin, driver.ode.t)
      write_data(sim)
    end
	if save_vtk_every > 0
		if i%save_vtk_every == 0
		  save_vtk(sim, @sprintf("%s_%d", sim.name, i))
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


function run_until(sim::AbstractSim, t_end::Float64; save_data=true)
      rk_data = sim.driver.ode
      if t_end < rk_data.t - rk_data.step
          println("Run_until: t_end >= rk_data.t - rk_data.step")
          return
      elseif t_end == rk_data.t
          rk_data.omega_t[:] = rk_data.omega[:]
          omega_to_spin(rk_data.omega_t, sim.prespin, sim.spin, sim.nxyz)
          sim.saver.t = t_end
          sim.saver.nsteps += 1
          if save_data
              compute_system_energy(sim, sim.spin, t_end)
              write_data(sim)
		  end
          return
      elseif t_end > rk_data.t - rk_data.step && rk_data.step > 0 && t_end < rk_data.t
          interpolation_dopri5(rk_data, t_end)
          omega_to_spin(rk_data.omega_t, sim.prespin, sim.spin, sim.nxyz)
          sim.saver.t = t_end
          sim.saver.nsteps += 1
          if save_data
              compute_system_energy(sim, sim.spin, t_end)
              write_data(sim)
          end
          return
      end

      # so we have t_end > self.t
      if rk_data.step_next<=0
          rk_data.step_next = compute_init_step(sim, t_end - rk_data.t)
      end

      while rk_data.t < t_end
          ratio = (t_end - rk_data.t)/rk_data.step_next
          if ratio<1.2 && ratio>0.8
              rk_data.step_next = t_end - rk_data.t
          end

          advance_step(sim, rk_data)
      end

      interpolation_dopri5(rk_data, t_end)
      omega_to_spin(rk_data.omega_t, sim.prespin, sim.spin, sim.nxyz)
      sim.saver.t = t_end
      sim.saver.nsteps += 1
	  if save_data
		  compute_system_energy(sim, sim.spin, t_end)
		  write_data(sim)
	  end
end
