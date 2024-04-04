using JLD2

export Sim, init_m0, set_Ms

"""
    Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="DormandPrince")

Create a simulation instance for given mesh.

"""
function Sim(mesh::FDMesh; driver="LLG", name="dyn", integrator="DormandPrince", save_data=true)

    sim = MicroSim()

    sim.name = name
    sim.mesh = mesh
    n_total = mesh.nx*mesh.ny*mesh.nz
    sim.n_total = n_total

    T = single_precision.x ? Float32 : Float64

    sim.spin = KernelAbstractions.zeros(backend[], T, 3*n_total) 
    sim.prespin = KernelAbstractions.zeros(backend[], T, 3*n_total) 
    sim.field = KernelAbstractions.zeros(backend[], T, 3*n_total) 
    sim.energy = KernelAbstractions.zeros(backend[], T, n_total) 
    sim.Ms = KernelAbstractions.zeros(backend[], T, n_total) 
    sim.pins = KernelAbstractions.zeros(backend[], Bool, n_total) 
    sim.driver_name = driver
    sim.driver = create_driver(driver, integrator, n_total)    
    sim.interactions = []
    sim.save_data = save_data
    sim.saver = create_saver(string(name, "_", lowercase(driver), ".txt"), driver)
   
   @info "Standard Sim (CPU) has been used."
   return sim
end

"""
    set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)

Set the saturation magnetization Ms of the studied system. For example,

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
    T = single_precision.x ? Float32 : Float64
    Ms_a = zeros(T, sim.n_total)
    init_scalar!(Ms_a, sim.mesh, Ms)

    if any(isnan, Ms_a)
      error("NaN is given by the input Ms!")
    end

    copyto!(sim.Ms, Ms_a)
    return true
end


"""
    set_Ms_cylindrical(sim::MicroSim, Ms::Number; axis=ez, r1=0, r2=0)

Set the saturation magnetization Ms of the studied system to a cylindrical shape,
axis could be JuMag.ex, JuMag.ey or JuMag.ez.
"""
function set_Ms_cylindrical(sim::MicroSim, Ms::Number; axis=ez, r1=0, r2=0)
    geo = create_cylinder(sim.mesh, axis, r1=r1, r2=r2)
    for i in 1:sim.n_total
        if geo.shape[i]
            sim.Ms[i] = Ms
        end
    end
    return true
end

"""
    set_Ms(sim::AbstractSim, geo::Geometry, Ms::Number)

Set the saturation magnetization Ms within the Geometry.
"""
function set_Ms(sim::AbstractSim, geo::Geometry, Ms::Number)
    update_scalar_geometry(sim.Ms, geo, Ms)
    return true
end


"""
    set_pinning(sim::MicroSim, ids::ArrayOrFunction)

Pinning the spins at the given ids.

```julia
function pinning_boundary(i,j,k,dx,dy,dz)
    if i == 1 || i == 100
        return true
    end
    return false
end
set_pinning(sim, pinning_boundary)
```
"""
function set_pinning(sim::MicroSim, ids::ArrayOrFunction)
    init_scalar!(sim.pins, sim.mesh, ids)
    return true
end

function set_ux(sim::AbstractSim, init_ux)
    init_scalar!(sim.driver.ux, sim.mesh, init_ux)
end

function set_ux_bounary(sim::AbstractSim, ux)
    return set_ux_bounary_implement(sim, ux)
end

function set_uy(sim::AbstractSim, init_uy)
    init_scalar!(sim.driver.uy, sim.mesh, init_uy)
end

function set_uz(sim::AbstractSim, init_uz)
    init_scalar!(sim.driver.uz, sim.mesh, init_uz)
end

function set_aj(sim::AbstractSim, init_aj)
	init_scalar!(sim.driver.aj, sim.mesh, init_aj)
end

function average_m(sim::AbstractSim)
  b = reshape(sim.spin, 3, sim.n_total)
  mx,my,mz = 0.0,0.0,0.0
  n = 0
  ms = isa(sim, MicroSim) ? sim.Ms : sim.mu_s
  for i = 1:sim.n_total
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

"""
    init_m0(sim::MicroSim, m0::TupleOrArrayOrFunction; norm=true)

Set the initial magnetization of the system. If `norm=false` the magnetization array will be not normalised.
Examples:

```julia
   init_m0(sim, (1,1,1))
```
or
```julia
   init_m0(sim, (1,1,1), norm=false)
```
or
```julia
   function uniform_m0(i,j,k,dx,dy,dz)
       return (0,0,1)
   end
   init_m0(sim, uniform_m0)
```
"""
function init_m0(sim::AbstractSim, m0::TupleOrArrayOrFunction; norm=true)

  spin = Array(sim.spin)

  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.n_total)
  end

  Ms = isa(sim.mesh, FDMesh) ? Array(sim.Ms) : Array(sim.mu_s)
  for i = 1:sim.n_total
      if Ms[i] == 0
          spin[3*i-2] = 0
          spin[3*i-1] = 0
          spin[3*i] = 0
      end
  end

  if any(isnan, spin)
    error("NaN is given by the input m0!")
  end

  copyto!(sim.spin, spin)
  copyto!(sim.prespin, sim.spin)
  return true
end


"""
    set_driver(sim::AbstractSim; driver="LLG", integrator="DormandPrince")

Set the driver of the simulation, can be used to switch the driver.
"""
function set_driver(sim::AbstractSim; driver="LLG", integrator="DormandPrince", args...)

    args = Dict(args)

    if sim.driver_name != driver
        # if the driver is updated, we create a new saver
        if sim.save_data
            sim.saver = create_saver(string(sim.name, "_", lowercase(driver), ".txt"), driver)
        end

        if isa(sim, AbstractSimGPU)
            sim.driver = create_driver_gpu(driver, integrator, sim.n_total)
        else
            sim.driver = create_driver(driver, integrator, sim.n_total)
        end
        sim.driver_name = driver    
    end

    # FIXME: we have to consider all the situations here
    
    if haskey(args, :alpha) && startswith(driver, "LLG")
        sim.driver.alpha = args[:alpha]
        delete!(args, :alpha)
    end

    if haskey(args, :gamma) && startswith(driver, "LLG")
        sim.driver.gamma = args[:gamma]
        delete!(args, :gamma)
    end

    if haskey(args, :beta) && startswith(driver, "LLG_STT")
        sim.driver.beta = args[:beta]
        delete!(args, :beta)
    end

    if haskey(args, :ux) && startswith(driver, "LLG_STT")
        set_ux(sim, args[:ux])
        delete!(args, :ux)
    end

    if haskey(args, :uy) && startswith(driver, "LLG_STT")
        set_uy(sim, args[:uy])
        delete!(args, :uy)
    end

    if haskey(args, :uz) && startswith(driver, "LLG_STT")
        set_uz(sim, args[:uz])
        delete!(args, :uz)
    end
end

"""
    relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, save_m_every = 10, save_ovf_every=-1, ovf_format = "binary", ovf_folder="ovfs", save_vtk_every=-1, vtk_folder="vtks",fields::Array{String, 1} = String[])

Relax the system using `LLG` or `SD` driver. This function works for both micromagnetic (FD and FE) and atomistic simulations, in both CPU and GPU. 

`maxsteps` is the maximum steps allowed to run. 

`stopping_dmdt` is the main stop condition, both for both for `LLG` and `SD` drivers. For standard micromagnetic simulaition, 
the typical value of `stopping_dmdt` is in the range of [0.01, 1].  In the `SD` driver, the time is not strictly defined. 
To make it comparable for the `LLG` driver, we multiply a factor of `gamma`. However, for the atomistic model 
with dimensionless unit, this factor should not be used. In this situation, `using_time_factor` should be set to `false`.

The magnetization (spins) can be stored in ovfs or vtks. ovf format can be chosen in "binary"(float64),"binary8"(float64), "binary4"(float32), "text"

Fields can be stored in vtks as well

```julia
relax(sim, save_vtk_every = 10, fields = ["demag", "exch", "anis"])
```
"""
function relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, using_time_factor=true, 
    save_m_every = -1, save_ovf_every=-1, save_vtk_every=-1, 
    ovf_format = "binary", ovf_folder="ovfs", vtk_folder="vtks", fields::Array{String, 1} = String[])

    # to dertermine which driver is used.
    llg_driver = false
    if isa(sim.driver, LLG) || (_cuda_available.x && isa(sim.driver, LLG_GPU))
        llg_driver = true
    end

    time_factor =  using_time_factor ? 2.21e5/2 : 1.0

    N_spins = sim.n_total

    if _cuda_available.x && (isa(sim, MicroSimGPU) || isa(sim, AtomicSimGPU))
        T = _cuda_using_double.x ? Float64 : Float32
        dm = CUDA.zeros(T, 3*sim.n_total)
    else
        dm = zeros(Float64,3*N_spins)
    end


    dmdt_factor = (2 * pi / 360) * 1e9
    if _cuda_available.x && isa(sim, AtomicSimGPU)
        dmdt_factor = 1.0
    end

    if save_ovf_every > 0
        isdir(ovf_folder) || mkdir(ovf_folder)
    end

    if save_vtk_every > 0
        isdir(vtk_folder) || mkdir(vtk_folder)
    end

    step = 0
    driver = sim.driver
    @info @sprintf("Running Driver : %s.", typeof(driver))
    for i=1:maxsteps

        run_step(sim, driver)

        step_size = llg_driver ? driver.ode.step : driver.tau/time_factor

        compute_dm!(dm, sim.prespin, sim.spin, N_spins)
        max_dmdt = maximum(dm)/step_size

        t = llg_driver ? sim.driver.ode.t : 0.0
        if llg_driver
            @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
                            i, step_size, t, max_dmdt/dmdt_factor)
        else
            @info @sprintf("step =%5d  step_size=%10.6e    max_dmdt=%10.6e",
                            i, step_size, max_dmdt/dmdt_factor)
        end
        
        if save_m_every>0 && i%save_m_every == 0
            compute_system_energy(sim, sim.spin, t)
            write_data(sim)
        end

        if save_vtk_every > 0 && i%save_vtk_every == 0
            save_vtk_points(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)), fields = fields)
        end

        if save_ovf_every > 0 && i%save_ovf_every == 0
            save_ovf(sim, joinpath(ovf_folder, @sprintf("%s_%d", sim.name, i)), dataformat = ovf_format)
        end
        
        if sim.save_data
            sim.saver.nsteps += 1
        end

        if max_dmdt < stopping_dmdt*dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
            if save_m_every>0 
                compute_system_energy(sim, sim.spin, t)
                write_data(sim)
            end
    
            if save_vtk_every > 0
                save_vtk_points(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)), fields = fields)
            end
    
            if save_ovf_every > 0 
                save_ovf(sim, joinpath(ovf_folder, @sprintf("%s_%d", sim.name, i)), dataformat = ovf_format)
            end
            step = i
            break
        end
    end

    if step == maxsteps
        if save_vtk_every > 0
            save_vtk_points(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)), fields = fields)
        end

        if save_ovf_every > 0 
            save_ovf(sim, joinpath(ovf_folder, @sprintf("%s_%d", sim.name, i)), dataformat = ovf_format)
        end
    end
  return nothing
end


function run_until(sim::AbstractSim, t_end::Float64, integrator::IntegratorCayley, save_data::Bool)
      if t_end < integrator.t - integrator.step
          println("Run_until: t_end >= integrator.t - integrator.step")
          return
      elseif t_end == integrator.t
          integrator.omega_t[:] = integrator.omega[:]
          omega_to_spin(integrator.omega_t, sim.prespin, sim.spin, sim.n_total)
          sim.saver.t = t_end
          sim.saver.nsteps += 1
          if save_data
              compute_system_energy(sim, sim.spin, t_end)
              write_data(sim)
        end
          return
      elseif t_end > integrator.t - integrator.step && integrator.step > 0 && t_end < integrator.t
          interpolation_dopri5(integrator, t_end)
          omega_to_spin(integrator.omega_t, sim.prespin, sim.spin, sim.n_total)
          sim.saver.t = t_end
          sim.saver.nsteps += 1
          if save_data
              compute_system_energy(sim, sim.spin, t_end)
              write_data(sim)
          end
          return
      end

      # so we have t_end > self.t
      if integrator.step_next<=0
          integrator.step_next = compute_init_step(sim, t_end - integrator.t)
      end

      while integrator.t < t_end
          ratio = (t_end - integrator.t)/integrator.step_next
          if ratio<1.2 && ratio>0.8
              integrator.step_next = t_end - integrator.t
          end

          advance_step(sim, integrator)
      end

      interpolation_dopri5(integrator, t_end)
      omega_to_spin(integrator.omega_t, sim.prespin, sim.spin, sim.n_total)
      sim.saver.t = t_end
      sim.saver.nsteps += 1
      if save_data
          compute_system_energy(sim, sim.spin, t_end)
          write_data(sim)
      end
      return nothing
end


function run_until(sim::AbstractSim, t_end::Float64, integrator::Integrator, save_data::Bool)
    if t_end < integrator.t - integrator.step
        @info("Run_until: t_end >= integrator.t - integrator.step")
        return
    elseif t_end == integrator.t
        sim.saver.t = t_end
        sim.saver.nsteps += 1
        if save_data
            compute_system_energy(sim, sim.spin, t_end)
            write_data(sim)
        end
        return
    end

    # so we have t_end > self.t
    if integrator.step_next<=0
        integrator.step_next = compute_init_step_DP(sim, t_end - integrator.t)
    end

    while integrator.t < t_end
        if integrator.step_next + integrator.t> t_end
            integrator.step_next = t_end - integrator.t
        end
        advance_step(sim, integrator)
    end

    sim.saver.t = t_end
    sim.saver.nsteps += 1
    if save_data
        compute_system_energy(sim, sim.spin, t_end)
        write_data(sim)
    end
end

function run_until(sim::AbstractSim, t_end::Float64; save_data=true)
    run_until(sim, t_end, sim.driver.ode, save_data)
end

"""
    create_sim(mesh; args...)

Create a micromagnetic simulation instance with given arguments. 

- `mesh`: a mesh has to be provided to start the simulation. The mesh could be [`FDMesh`](@ref), 
[`FDMeshGPU`](@ref), [`CubicMeshGPU`](@ref), or [`TriangularMeshGPU`](@ref).

# Arguments
- `name` : the simulation name, should be a string.
- `driver` : the driver name, should be a string. By default, the driver is "SD".
- `alpha` : the Gilbert damping in the LLG equation, should be a number.
- `beta` : the nonadiabatic strength in the LLG equation with spin transfer torques (zhang-li model), should be a number.
- `gamma` : the gyromagnetic ratio, default value = 2.21e5.
- `ux`, `uy` or `uz`: the strengths of the spin transfer torque.
- `Ms`: the saturation magnetization, should be [`NumberOrArrayOrFunction`](@ref). By default, Ms=8e5
- `mu_s`: the magnetic moment, should be [`NumberOrArrayOrFunction`](@ref). By default, mu_s=2*mu_B
- `A` or `J`: the exchange constant, should be [`NumberOrArrayOrFunction`](@ref).
- `D` : the DMI constant, should be [`NumberOrArrayOrFunction`](@ref).
- `dmi_type` : the type of DMI, could be "bulk" or "interfacial".
- `Ku`: the anisotropy constant, should be [`NumberOrArrayOrFunction`](@ref).
- `axis`: the anisotropy axis, should be a tuple, such as (0,0, 1)
- `demag` : include demagnetization or not, should be a boolean, i.e., true or false.
- `H`: the external field, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref). 
- `m0` : the initial magnetization, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref). 
"""
function create_sim(mesh; args...)
    #convert args to a dict
    args = Dict(args)
    #Ms=8e5, A=1.3e-11, D=0, Ku=0, axis=(0,0,1), H=(0,0,0), m0=(0,0,1), demag=false, driver="SD", name="relax"
    driver = haskey(args, :driver) ? args[:driver] : "SD"
    name = haskey(args, :name) ? args[:name] : "unnamed"

    #Create the mesh using given driver and name
    sim = Sim(mesh, driver=driver, name = name)

    #If the simulation is the standard micromagnetic simulation.
    if isa(mesh, FDMesh) || isa(mesh, FDMeshGPU) || isa(mesh, FDMeshGPU)

        # we set the Ms anyway
        Ms = haskey(args, :Ms) ? args[:Ms] : 8e5
        set_Ms(sim, Ms)

        # add the exchange if A is given
        if haskey(args, :A)
            add_exch(sim, args[:A])
        end

        # add the DMI if D is given
        if haskey(args, :D)
            dmi_type = haskey(args, :dmi_type) ? args[:dmi_type] : "bulk"
            add_dmi(sim, args[:D], type=dmi_type)
        end

        # add the demag
        if haskey(args, :demag) && args[:demag]
            add_demag(sim)
        end

        for key in [:Ms, :A, :D, :demag]
            haskey(args, key) && delete!(args, key)
        end

    
    #If the simulation is atomistic
    elseif isa(mesh, AtomicMeshGPU)

        mu_s = haskey(args, :mu_s) ? args[:mu_s] : 2*mu_B
        set_mu_s(sim, mu_s)

        # add the exchange if A is given
        if haskey(args, :J)
            add_exch(sim, args[:J])
        end
        
        # add the DMI if D is given
        if haskey(args, :D)
            add_dmi(sim, args[:D])
        end

        for key in [:mu_s, :J, :D]
            haskey(args, key) && delete!(args, key)
        end
    else
        error("This info is for debug.")
    end

    # add the anisotropy
    if haskey(args, :Ku)
        axis = haskey(args, :axis) ? args[:axis] : (0,0,1)
        add_anis(sim, args[:Ku], axis=axis)
        haskey(args, :axis) && delete!(args, :axis)
    end

    # add the external field
    if haskey(args, :H)
        add_zeeman(sim, args[:H])
    end

    if haskey(args, :alpha) && startswith(driver, "LLG")
        sim.driver.alpha = args[:alpha]
        delete!(args, :alpha)
    end

    if haskey(args, :gamma) && startswith(driver, "LLG")
        sim.driver.gamma = args[:gamma]
        delete!(args, :gamma)
    end

    if haskey(args, :beta) && startswith(driver, "LLG_STT")
        sim.driver.beta = args[:beta]
        delete!(args, :beta)
    end

    if haskey(args, :ux) && startswith(driver, "LLG_STT")
        set_ux(sim, args[:ux])
        delete!(args, :ux)
    end

    if haskey(args, :uy) && startswith(driver, "LLG_STT")
        set_uy(sim, args[:uy])
        delete!(args, :uy)
    end

    if haskey(args, :uz) && startswith(driver, "LLG_STT")
        set_uz(sim, args[:uz])
        delete!(args, :uz)
    end

    # set m0 anyway
    m0_value = haskey(args, :m0) ? args[:m0] : (0.8, 0.6, 0)
    init_m0(sim, m0_value)
  
    for key in [:driver, :name, :Ku, :H, :m0]
        haskey(args, key) && delete!(args, key)
    end
    
    for key in args
        @warn @sprintf("Key '%s' is not used.", key)
    end
    
    return sim
    
  end


"""
    run_sim(sim::AbstractSim; steps=10, dt=1e-12, save_data=true, save_m_every=1, saver=nothing)

Run the simulation to the time `steps*dt`.

- steps : the total steps of the simulation
- dt : the time interval of each step, so the total simulation time is `steps*dt`
- save_data : saving the overall data such as energies and average magnetization of the simulation at each step
- save_m_every : save magnetization for every `save_m_every` step, a negative save_m_every will disable the magnetization saving.
- saver : a saver struct, by default it will use sim's saver. But you can use customized saver instead. For example, if we want to compute the guiding
center and save it to a text file, we can define the following saver
```julia
    customized_saver = init_saver("output.txt", "LLG")
    push!(custom_saver.items, SaverItem("center", "m", compute_guiding_center))
    run_sim(sim, saver=customized_saver)
```

"""
function run_sim(sim::AbstractSim; steps=10, dt=1e-12, save_data=true, save_m_every=1, saver=nothing)

    if save_data && saver===nothing
        saver = sim.saver
    end

    file = jldopen(@sprintf("%s.jdl2", sim.name), "w")

    file["mesh/nx"] = sim.mesh.nx
    file["mesh/ny"] = sim.mesh.ny
    file["mesh/nz"] = sim.mesh.nz
    file["mesh/dx"] = sim.mesh.dx
    file["mesh/dy"] = sim.mesh.dy
    file["mesh/dz"] = sim.mesh.dz

    file["steps"] = steps
    file["dt"] = dt
    file["save_m_every"] = save_m_every

    if save_m_every > 0
        m_group = JLD2.Group(file, "m")
    end

    for i = 1:steps
        run_until(sim, i*dt, save_data=save_data)
        @info @sprintf("step =%5d  t = %10.6e", i, i*dt)

        save_data && write_data(sim, saver)

        if (save_m_every>1 && i%save_m_every == 1) || save_m_every == 1
            index = @sprintf("%d", i)
            m_group[index] = Array(sim.spin)
        end
    end
    close(file)
end