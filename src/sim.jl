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

"""
    set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)

Set the saturation magnetization of the studied system. As can be guessed from the type union `NumberOrArrayOrFunction`,
Ms could be a number or an array or a function. If Ms is an array, its length should be equal to
the size of the system. If Ms is a function, it should take six parameters in the form
`(i,j,k,dx,dy,dz)` where `i,j,k` are the indices of the mesh cells and `dx,dy,dz` are the cellsizes.
For example,

```julia
   set_Ms(sim, 8.6e5)
```
or
```julia
function circular_Ms(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8.6e5
    end
    return 0.0
end
set_Ms(sim, circular_Ms)
```

"""
function set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)
    init_scalar!(sim.Ms, sim.mesh, Ms)
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

function init_m0(sim::AtomicSim, m0::TupleOrArrayOrFunction; norm=true)
  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
    normalise(sim.prespin, sim.nxyz)
  end
  for i = 1:sim.nxyz
      if sim.mu_s[i] == 0.0
          sim.prespin[3*i-2] = 0
          sim.prespin[3*i-1] = 0
          sim.prespin[3*i] = 0
      end
  end
  sim.spin[:] .= sim.prespin[:]
end

"""
    init_m0(sim::MicroSim, m0::TupleOrArrayOrFunction; norm=true)

Set the initial magnetization of the system, the input `m0` could be a Tuple with length 3
or an array with length 3*N where N is the system size or a function with six parameters `(i,j,k,dx,dy,dz)`.
"""
function init_m0(sim::MicroSim, m0::TupleOrArrayOrFunction; norm=true)
  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
    normalise(sim.prespin, sim.nxyz)
  end
  for i = 1:sim.nxyz
      if sim.Ms[i] == 0.0
          sim.prespin[3*i-2] = 0
          sim.prespin[3*i-1] = 0
          sim.prespin[3*i] = 0
      end
  end
  sim.spin[:] .= sim.prespin[:]
end

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Add a static Zeeman energy to the simulation.
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
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

"""
    update_zeeman(z::Zeeman, H0::Tuple)

Set the external field of Zeeman z to H0 where H0 is a tuple.
For example,
```julia
   ze = add_zeeman(sim, (0,0,0)) #create a zeeman energy with field (0,0,0) A/m
   update_zeeman(ze, (0,0,1e5)) #change the field to (0,0,1e5) A/m
```
"""
function update_zeeman(z::Zeeman, H0::Tuple)
  b = reshape(z.field, 3, Int(length(z.field)/3))
  b[1, :] .= H0[1]
  b[2, :] .= H0[2]
  b[3, :] .= H0[3]
  return z
end

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")

Add a time varying zeeman to system. The input `H0` could be

 - A tuple with length 3, for example, `(0,0,1e2)`.
 - An array with length `3N` where N is the number of total spins.
 - A function with six parameters `(i,j,k,dx,dy,dz)` and return a tuple, for example,
   ```julia
    function linear_Hz(i, j, k, dx, dy, dz)
       nx = 100
       Hz = i/nx*1000
       return (0, 0, Hz)
    end
   ```

The input `ft` is a function of time `t` and its return value should be a tuple. For example,

```julia
  function time_fun(t)
    return (sin(t), cos(t), 0)
  end
```
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")
  nxyz = sim.nxyz
  init_field = zeros(Float64, 3*nxyz)
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  init_vector!(init_field, sim.mesh, H0)

  zeeman =  TimeZeeman(ft, init_field, field, energy, name)
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


function add_exch(sim::AbstractSim, A::NumberOrArrayOrFunction; name="exch")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  Spatial_A = zeros(Float64, nxyz)
  init_scalar!(Spatial_A , sim.mesh, A)
  if isa(sim, MicroSim)
    exch = Exchange(Spatial_A , field, energy, name)
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


function add_exch_rkky(sim::AbstractSim, sigma::Float64, Delta::Float64; name="rkky")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  exch = ExchangeRKKY(sigma, Delta, field, energy, name)

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

"""
    add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)

Add Demag to the system. ``N_x``, ``N_y`` and ``N_z`` can be used to describe the macro boundary conditions which means that
the given mesh is repeated ``2N_x+1``, ``2N_y+1`` and ``2N_z+1`` times in ``x``, ``y`` and ``z`` direction, respectively.
"""
function add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)
  demag = init_demag(sim, Nx, Ny, Nz)
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
  lt=sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
  naxis=(axis[1]^2/lt,axis[2]^2/lt,axis[3]^2/lt)
  anis =  Anisotropy(Kus, naxis, field, energy, name)
  push!(sim.interactions, anis)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  return anis
end

function relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, stopping_torque=0.1, save_m_every = 10, save_vtk_every=-1, vtk_folder="vtks")
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
      relax_llg(sim, maxsteps, Float64(stopping_dmdt), save_m_every, save_vtk_every, vtk_folder)
  else
      relax_energy(sim, maxsteps, Float64(stopping_torque), save_m_every, save_vtk_every, vtk_folder)
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
  for i=1:maxsteps
    run_step(sim, driver)
    abs!(gk_abs, driver.gk)  #max_torque = maximum(abs.(driver.gk)) eats gpu memory???
    max_torque = maximum(gk_abs)
	@info @sprintf("step=%5d  tau=%10.6e  max_torque=%10.6e", i, driver.tau, max_torque)
    if save_m_every>0 && i%save_m_every == 0
      compute_system_energy(sim, sim.spin, 0.0)
      write_data(sim)
    end
    if save_vtk_every > 0 && i%save_vtk_every == 0
        save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
    end
    sim.saver.nsteps += 1
    if max_torque < stopping_torque
      @info @sprintf("max_torque (mxmxH) is less than stopping_torque=%g, Done!", stopping_torque)
      if save_m_every>0
          compute_system_energy(sim, sim.spin, 0.0)
          write_data(sim)
      end
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
  for i=1:maxsteps
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
    if save_m_every>0 && i%save_m_every == 0
      compute_system_energy(sim, sim.spin, sim.driver.ode.t)
      write_data(sim)
    end
    if save_vtk_every > 0 && i%save_vtk_every == 0
        save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
    end
    sim.saver.t = rk_data.t
    sim.saver.nsteps += 1
    if max_dmdt < stopping_dmdt*dmdt_factor
      @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
      if save_m_every>0
          compute_system_energy(sim, sim.spin, sim.driver.ode.t)
          write_data(sim)
      end
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
