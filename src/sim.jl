"""
    Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="Dopri5")

Create a simulation instance for given mesh.
"""
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

function set_aj(sim::AbstractSim, init_aj)
	init_scalar!(sim.driver.aj, sim.mesh, init_aj)
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

  zeeman =  Zeeman(field, energy, name)
  push!(sim.interactions, zeeman)

  if isa(H0, Tuple)
      push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
      push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
      id = length(sim.interactions)
      fun = o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
      push!(sim.saver.results, fun)
  end

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  return zeeman
end

function add_zeeman(sim::AbstractSim, H0::Any, funs::Tuple{Function,Function,Function}; name="timezeeman")
  nxyz = sim.nxyz
  init_field = zeros(Float64, 3*nxyz)
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  init_vector!(init_field, sim.mesh, H0)

  zeeman =  TimeZeeman(funs[1], funs[2], funs[3], init_field, field, energy, name)
  push!(sim.interactions, zeeman)

  if isa(H0, Tuple)
      push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
      push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
      id = length(sim.interactions)
      fun = o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
      push!(sim.saver.results, fun)
  end

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

function relax(sim::AbstractSim; maxsteps=10000, init_step = 1e-13, stopping_dmdt=0.01, stopping_torque=0.1, save_m_every = 10, save_vtk_every=-1, vtk_folder="vtks")
  is_relax_llg = false
  if _using_gpu.x && isa(sim.driver, LLG_GPU)
        is_relax_llg = true
  end

  if isa(sim.driver, LLG)
      is_relax_llg = true
  end

  if !isdir(vtk_folder)
      mkdir(vtk_folder)
  end

  if is_relax_llg
      relax_llg(sim, maxsteps, stopping_dmdt, save_m_every, save_vtk_every, vtk_folder)
  else
      relax_energy(sim, maxsteps, stopping_torque, save_m_every, save_vtk_every, vtk_folder)
  end
  return nothing
end


function relax_energy(sim::AbstractSim, maxsteps::Int64, stopping_torque::Float64,
                     save_m_every::Int64, save_vtk_every::Int64, vtk_folder::String)

  if _using_gpu.x && isa(sim, MicroSimGPU)
      T = _cuda_using_double.x ? Float64 : Float32
      gk_abs = cuzeros(T, 3*sim.nxyz)
  else
      gk_abs = zeros(Float64,3*sim.nxyz)
  end

  driver = sim.driver
  for i=0:maxsteps-1
    run_step(sim, driver)
    abs!(gk_abs, driver.gk)  #max_torque = maximum(abs.(driver.gk)) eats gpu memory???
    max_torque = maximum(gk_abs)
	@info @sprintf("step=%5d  tau=%10.6e  max_torque=%10.6e", i, driver.tau, max_torque)
    if i%save_m_every == 0
      compute_system_energy(sim, sim.spin, 0.0)
      write_data(sim)
    end
    if save_vtk_every > 0
        if i%save_vtk_every == 0
            save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
        end
    end
    sim.saver.nsteps += 1
    if max_torque < stopping_torque
      @info @sprintf("max_torque (mxmxH) is less than stopping_torque=%g, Done!", stopping_torque)
      compute_system_energy(sim, sim.spin, 0.0)
      write_data(sim)
      if save_vtk_every > 0
          save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
      end
      break
    end
  end
  return nothing
end

function relax_llg(sim::AbstractSim, maxsteps::Int64, stopping_dmdt::Float64,
         save_m_every::Int64, save_vtk_every::Int64, vtk_folder::String)
  step = 0
  rk_data = sim.driver.ode

  dmdt_factor = 1.0
  if isa(sim, MicroSim) || (_using_gpu.x && isa(sim, MicroSimGPU))
    dmdt_factor = (2 * pi / 360) * 1e9
  end
  for i=0:maxsteps-1
    advance_step(sim, rk_data)
    step_size = rk_data.step
    #omega_to_spin(rk_data.omega, sim.prespin, sim.spin, sim.nxyz)
    max_dmdt = 0.0
    if _using_gpu.x
        compute_dm(rk_data.omega_t, sim.prespin, sim.spin, sim.nxyz)
        max_dmdt = maximum(rk_data.omega_t)/step_size
    else
        max_dmdt = compute_dmdt(sim.prespin, sim.spin, sim.nxyz, step_size)
    end

    @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
                   i, rk_data.step, rk_data.t, max_dmdt/dmdt_factor)
    if i%save_m_every == 0
      compute_system_energy(sim, sim.spin, sim.driver.ode.t)
      write_data(sim)
    end
    if save_vtk_every > 0
        if i%save_vtk_every == 0
            save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
        end
    end
    sim.saver.t = rk_data.t
    sim.saver.nsteps += 1
    if max_dmdt < stopping_dmdt*dmdt_factor
      @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
      compute_system_energy(sim, sim.spin, sim.driver.ode.t)
      write_data(sim)
      if save_vtk_every > 0
          save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
      end
      break
    end
  end
  return nothing
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
      return nothing
end
