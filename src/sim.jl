using JLD2

export Sim, init_m0, set_Ms, set_Ms_cylindrical, run_until, relax, create_sim, run_sim,
       set_driver, set_pinning, set_ux, set_uy, set_uz

"""
    Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="DormandPrince")

Create a simulation instance for given mesh.

"""
function Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="DormandPrince",
             save_data=true)
    T = Float[]

    if isa(mesh, FDMesh)
        sim = MicroSim{T}()
    else
        sim = AtomisticSim{T}()
    end

    sim.name = name
    sim.mesh = mesh
    n_total = mesh.nx * mesh.ny * mesh.nz
    sim.n_total = n_total
    sim.spin = create_zeros(3 * n_total)
    sim.prespin = create_zeros(3 * n_total)
    sim.field = create_zeros(3 * n_total)
    sim.energy = create_zeros(n_total)
    if isa(mesh, FDMesh)
        sim.mu0_Ms = create_zeros(n_total)
    else
        sim.mu_s = create_zeros(n_total)
    end
    sim.pins = create_zeros(Bool, n_total)
    sim.driver_name = driver
    sim.driver = create_driver(driver, integrator, n_total)
    sim.interactions = []
    sim.save_data = save_data
    sim.saver = create_saver(string(name, "_", lowercase(driver), ".txt"), driver)

    if isa(mesh, FDMesh)
        @info "MicroSim has been created."
    else
        @info "AtomisticSim has been created."
    end

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
    T = Float[]
    Ms_a = zeros(T, sim.n_total)
    init_scalar!(Ms_a, sim.mesh, Ms)

    if any(isnan, Ms_a)
        error("NaN is given by the input Ms!")
    end

    Ms_a .*= mu_0  #we convert A/m to Tesla

    copyto!(sim.mu0_Ms, Ms_a)
    return true
end

"""
    set_Ms(sim::AbstractSim, geo::Shape, Ms::Number)

Set the saturation magnetization Ms within the Shape.
"""
function set_Ms(sim::AbstractSim, shape::Shape, Ms::Number)
    init_scalar!(sim.mu0_Ms, sim.mesh, shape, Ms * mu_0)
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
    pins = zeros(Bool, sim.n_total)
    init_scalar!(pins, sim.mesh, ids)
    copyto!(sim.pins, pins)
    return true
end

function set_ux(sim::AbstractSim, init_ux)
    return init_scalar!(sim.driver.ux, sim.mesh, init_ux)
end

function set_ux_bounary(sim::AbstractSim, ux)
    return set_ux_bounary_implement(sim, ux)
end

function set_uy(sim::AbstractSim, init_uy)
    return init_scalar!(sim.driver.uy, sim.mesh, init_uy)
end

function set_uz(sim::AbstractSim, init_uz)
    return init_scalar!(sim.driver.uz, sim.mesh, init_uz)
end

function set_aj(sim::AbstractSim, init_aj)
    return init_scalar!(sim.driver.aj, sim.mesh, init_aj)
end

"""
Compute the average magnetization defined as

TODO: add the equations
"""
function average_m(sim::AbstractSim)
    b = reshape(sim.spin, 3, sim.n_total)
    ms = isa(sim, MicroSim) ? sim.mu0_Ms : sim.mu_s
    return Tuple(sum(b .* ms'; dims=2) ./ sum(ms))
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

    Ms = isa(sim.mesh, FDMesh) ? Array(sim.mu0_Ms) : Array(sim.mu_s)
    for i in 1:(sim.n_total)
        if Ms[i] == 0
            spin[3 * i - 2] = 0
            spin[3 * i - 1] = 0
            spin[3 * i] = 0
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
    set_driver(sim::AbstractSim, args::Dict)

Set the driver of the simulation. This function is not intended for users but for developers.
"""
function set_driver(sim::AbstractSim, args::Dict)
    
    # FIXME: we have to consider all the situations here

    # driver = string(typeof(sim.driver))
    driver = sim.driver_name

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

    if haskey(args, :ufun) && startswith(driver, "LLG_STT")
        sim.driver.ufun = args[:ufun]
        delete!(args, :ufun)
    end
end

"""
    relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, save_data_every=1, save_m_every=-1, using_time_factor=true)

Relax the system using `LLG` or `SD` driver. This function works for both micromagnetic and atomistic simulations.

`maxsteps` is the maximum steps allowed to run. 

`stopping_dmdt` is the main stop condition, both for both for `LLG` and `SD` drivers. For standard micromagnetic simulaition, 
the typical value of `stopping_dmdt` is in the range of [0.01, 1].  In the `SD` driver, the time is not strictly defined. 
To make it comparable for the `LLG` driver, we multiply a factor of `gamma`. However, for the atomistic model 
with dimensionless unit, this factor should not be used. In this situation, `using_time_factor` should be set to `false`.

`save_data_every` set the step for overall data saving such as energies and average magnetization. A negative `save_data_every` will
disable the data saving (`save_data_every=-1` will enable the data saving at the end of relaxing).

`save_m_every` set the step for magnetization saving, a negative `save_m_every` will disable the magnetization saving.

Examples:

```julia
    relax(sim, maxsteps=10000, stopping_dmdt=0.1)
```
"""
function relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, save_data_every=1, save_m_every=-1, using_time_factor=true)

    # to dertermine which driver is used.
    llg_driver = isa(sim.driver, LLG)

    time_factor = using_time_factor ? 2.21e5 / 2 : 1.0
    dmdt_factor = using_time_factor ? (2 * pi / 360) * 1e9 : 1

    file = jldopen(@sprintf("%s.jld2", sim.name), "w")
    file["mesh/nx"] = sim.mesh.nx
    file["mesh/ny"] = sim.mesh.ny
    file["mesh/nz"] = sim.mesh.nz
    file["mesh/dx"] = sim.mesh.dx
    file["mesh/dy"] = sim.mesh.dy
    file["mesh/dz"] = sim.mesh.dz

    file["steps"] = maxsteps
    file["save_m_every"] = save_m_every

    if save_m_every > 0
        m_group = JLD2.Group(file, "m")
    end

    N_spins = sim.n_total
    dm = create_zeros(3 * N_spins)

    driver = sim.driver
    @info @sprintf("Running Driver : %s.", typeof(driver))
    for i in 0:maxsteps
        @timeit timer "run_step" run_step(sim, driver)

        step_size = llg_driver ? driver.integrator.step : driver.tau / time_factor

        compute_dm!(dm, sim.prespin, sim.spin, N_spins)
        max_dmdt = maximum(dm) / step_size

        t = llg_driver ? sim.driver.integrator.t : 0.0
        if llg_driver
            @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e", i,
                           step_size, t, max_dmdt / dmdt_factor)
        else
            @info @sprintf("step =%5d  step_size=%10.6e  max_dmdt=%10.6e", i, step_size,
                           max_dmdt / dmdt_factor)
        end

        if save_data_every > 0 && i % save_data_every == 0
            compute_system_energy(sim, sim.spin, t)
            write_data(sim)
        end

        if sim.save_data
            sim.saver.nsteps += 1
        end

        if save_m_every > 0 && i % save_m_every == 0
            index = @sprintf("%d", i)
            m_group[index] = Array(sim.spin)
        end

        if max_dmdt < stopping_dmdt * dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
            if save_data_every > 0 || save_data_every == -1
                compute_system_energy(sim, sim.spin, t)
                write_data(sim)
            end

            if save_m_every > 0
                index = @sprintf("%d", i)
                m_group[index] = Array(sim.spin)
            end
            break
        end
    end

    close(file)

    return nothing
end

function run_until(sim::AbstractSim, t_end::Float64, integrator::IntegratorCayley,
                   save_data::Bool)
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
    elseif t_end > integrator.t - integrator.step &&
           integrator.step > 0 &&
           t_end < integrator.t
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
    if integrator.step_next <= 0
        integrator.step_next = compute_init_step(sim, t_end - integrator.t)
    end

    while integrator.t < t_end
        ratio = (t_end - integrator.t) / integrator.step_next
        if ratio < 1.2 && ratio > 0.8
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

function run_until(sim::AbstractSim, t_end::Float64, integrator::Integrator,
                   save_data::Bool)
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
    if integrator.step_next <= 0
        integrator.step_next = compute_init_step_DP(sim, t_end - integrator.t)
    end

    while integrator.t < t_end
        if integrator.step_next + integrator.t > t_end
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
    return nothing
end

function run_until(sim::AbstractSim, t_end::Float64; save_data=true)
    @timeit timer "run_until" run_until(sim, t_end, sim.driver.integrator, save_data)
    return nothing
end

"""
    set_driver(sim::AbstractSim; driver="LLG", integrator="DormandPrince", args...)

Set the driver of the simulation, can be used to switch the driver.
"""
function set_driver(sim::AbstractSim; driver="LLG", integrator="DormandPrince", args...)
    args = Dict(args)

    @info(@sprintf("The driver %s is used!", driver))

    # FIXME : Does not work if only the integrator changes  
    if sim.driver_name != driver
        # if the driver is updated, we create a new saver
        if sim.save_data
            sim.saver = create_saver(string(sim.name, "_", lowercase(driver), ".txt"),
                                     driver)
        end

        sim.driver = create_driver(driver, integrator, sim.n_total)
        sim.driver_name = driver
    end

    set_driver(sim, args)
end

"""
    create_sim(mesh; args...)

Create a micromagnetic simulation instance with given arguments. 

- `mesh`: a mesh has to be provided to start the simulation. The mesh could be [`FDMesh`](@ref), [`CubicMesh`](@ref), or [`TriangularMesh`](@ref).

# Arguments
- `name` : the simulation name, should be a string.
- `driver` : the driver name, should be a string. By default, the driver is "SD".
- `alpha` : the Gilbert damping in the LLG equation, should be a number.
- `beta` : the nonadiabatic strength in the LLG equation with spin transfer torques (zhang-li model), should be a number.
- `gamma` : the gyromagnetic ratio, default value = 2.21e5.
- `ux`, `uy` or `uz`: the strengths of the spin transfer torque.
- `ufun` : the time-dependent function for `u`. 
- `Ms`: the saturation magnetization, should be [`NumberOrArrayOrFunction`](@ref). By default, Ms=8e5
- `mu_s`: the magnetic moment, should be [`NumberOrArrayOrFunction`](@ref). By default, mu_s=2*mu_B
- `A` or `J`: the exchange constant, should be [`NumberOrArrayOrFunction`](@ref).
- `D` : the DMI constant, should be [`NumberOrArrayOrFunction`](@ref).
- `dmi_type` : the type of DMI, could be "bulk" or "interfacial".
- `Ku`: the anisotropy constant, should be [`NumberOrArrayOrFunction`](@ref).
- `axis`: the anisotropy axis, should be a tuple, such as (0,0, 1)
- `demag` : include demagnetization or not, should be a boolean, i.e., true or false. By default,  demag=false.
- `H`: the external field, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref). 
- `m0` : the initial magnetization, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref). 
- `T` : the temperature, should be should be [`NumberOrArrayOrFunction`](@ref).
- `shape` : the shape defines the geometry of the sample, where parameters are configured.
"""
function create_sim(mesh; args...)
    #convert args to a dict
    args = Dict(args)
    #Ms=8e5, A=1.3e-11, D=0, Ku=0, axis=(0,0,1), H=(0,0,0), m0=(0,0,1), demag=false, driver="SD", name="relax"
    driver = haskey(args, :driver) ? args[:driver] : "SD"
    name = haskey(args, :name) ? args[:name] : "unnamed"
    shape = haskey(args, :shape) ? args[:shape] : nothing

    #Create the mesh using given driver and name
    sim = Sim(mesh; driver=driver, name=name)

    #If the simulation is the standard micromagnetic simulation.
    if isa(mesh, FDMesh)

        # we set the Ms anyway
        Ms = haskey(args, :Ms) ? args[:Ms] : 8e5
        if shape === nothing
            set_Ms(sim, Ms)
        else
            set_Ms(sim, shape, Ms)
        end

        # add the exchange if A is given
        if haskey(args, :A)
            add_exch(sim, args[:A])
        end

        # add the DMI if D is given
        if haskey(args, :D)
            dmi_type = haskey(args, :dmi_type) ? args[:dmi_type] : "bulk"
            add_dmi(sim, args[:D]; type=dmi_type)
        end

        # add the demag
        if haskey(args, :demag) && args[:demag]
            add_demag(sim)
        end

        for key in [:Ms, :A, :D, :demag]
            haskey(args, key) && delete!(args, key)
        end

        #If the simulation is atomistic
    elseif isa(mesh, AtomisticMesh)
        mu_s = haskey(args, :mu_s) ? args[:mu_s] : 2 * mu_B
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
        axis = haskey(args, :axis) ? args[:axis] : (0, 0, 1)
        add_anis(sim, args[:Ku]; axis=axis)
        haskey(args, :axis) && delete!(args, :axis)
    end

    # add the external field
    if haskey(args, :H)
        add_zeeman(sim, args[:H])
    end

    # add the thermal noise
    if haskey(args, :T)
        add_thermal_noise(sim, args[:T])
    end

    # set the driver with args
    set_driver(sim, args)

    # set m0 anyway
    m0_value = haskey(args, :m0) ? args[:m0] : (0.8, 0.6, 0)
    init_m0(sim, m0_value)

    for key in [:driver, :name, :Ku, :H, :m0, :shape, :T]
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

- `steps` : the total steps of the simulation
- `dt`` : the time interval of each step, so the total simulation time is `steps*dt`
- `save_data` : saving the overall data such as energies and average magnetization of the simulation at each step
- `save_m_every` : save magnetization for every `save_m_every` step, a negative save_m_every will disable the magnetization saving.
- `saver` : a saver struct, by default it will use sim's saver. But you can use customized saver instead. For example, if we want to compute the guiding
center and save it to a text file, we can define the following saver
```julia
    customized_saver = init_saver("output.txt", "LLG")
    push!(custom_saver.items, SaverItem("center", "m", compute_guiding_center))
    run_sim(sim, saver=customized_saver)
```

"""
function run_sim(sim::AbstractSim; steps=10, dt=1e-12, save_data=true, save_m_every=1,
                 saver=nothing, call_back=nothing)
    if save_data && saver === nothing
        saver = sim.saver
    end

    if save_m_every > 0
        file = jldopen(@sprintf("%s.jld2", sim.name), "w")

        file["mesh/nx"] = sim.mesh.nx
        file["mesh/ny"] = sim.mesh.ny
        file["mesh/nz"] = sim.mesh.nz
        file["mesh/dx"] = sim.mesh.dx
        file["mesh/dy"] = sim.mesh.dy
        file["mesh/dz"] = sim.mesh.dz
    
        file["steps"] = steps
        file["dt"] = dt
        file["save_m_every"] = save_m_every

        m_group = JLD2.Group(file, "m")
    end

    for i in 0:steps
        run_until(sim, i * dt; save_data=save_data)
        @info @sprintf("step =%5d  t = %10.6e", i, i * dt)

        save_data && write_data(sim, saver)

        call_back !== nothing && call_back(sim, i * dt)

        if save_m_every > 0 && i % save_m_every == 0
            index = @sprintf("%d", i)
            m_group[index] = Array(sim.spin)
        end
    end

    if save_m_every > 0
        close(file)
    end
    return nothing
end
