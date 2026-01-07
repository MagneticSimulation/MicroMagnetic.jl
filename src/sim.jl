using KernelAbstractions
using Printf

export Sim, init_m0, set_Ms, set_Ms_cylindrical, run_until, relax, create_sim, run_sim,
       set_driver, set_pinning, set_ux, set_uy, set_uz, sim_with, advance_step, set_alpha

"""
    Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="DormandPrince",
             save_data=true)

Create a simulation instance for the given mesh with specified driver and integrator.

# Arguments
- `mesh::Mesh`: The computational mesh defining the simulation geometry. Can be:
  - `FDMesh`: Finite difference mesh for micromagnetics
  - `FEMesh`: Finite element mesh for micromagnetics  
  - `AtomisticMesh`: Creates atomistic simulation

# Keyword Arguments
- `driver::String="LLG"`: The simulation driver type. Options:
  - `"SD"`: Energy minimization (Steepest Descent)
  - `"LLG"`: Landau-Lifshitz-Gilbert equation
  - `"InertialLLG"` : Inertial LLG Equation
  - `"SpatialLLG"` : Spatial LLG equation allowing spatial damping constant.
  - `"LLG_STT"`: LLG with spin transfer torque
  - `"LLG_CPP"`: LLG with CPP spin transfer torque
- `name::String="dyn"`: Name identifier for the simulation
- `integrator::String="DormandPrince"`: Time integration method. Options:
  - Fixed step methods: `"Heun"`, `"RungeKutta"`, `"RungeKuttaCayley"`
  - Adaptive step methods: `"DormandPrince"` (DOPRI54), `"BS23"`, `"CashKarp54"`, `"Fehlberg54"` (RKF54)
  - Cayley-transform methods: `"DormandPrinceCayley"`, `"RungeKuttaCayley"`
- `save_data::Bool=true`: Whether to enable data saving during simulation

# Returns
- `MicroSim{Float}`: For FDMesh simulations
- `MicroSimFE{Float}`: For FEMesh simulations  
- `AtomisticSim{Float}`: For AtomisticMesh types

# Examples
```julia
# Create LLG simulation with Dormand-Prince integrator
mesh = FDMesh(nx=100, ny=100, nz=1)
sim = Sim(mesh, driver="LLG", integrator="DormandPrince")

# Create energy minimization simulation  
sim = Sim(mesh, driver="SD")
```
"""
function Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="DormandPrince",
             save_data=true)
    T = Float[]

    if isa(mesh, FDMesh)
        sim = MicroSim{T}()
    elseif isa(mesh, FEMesh)
        sim = MicroSimFE{T}()
    else
        sim = AtomisticSim{T}()
    end

    sim.name = name
    sim.mesh = mesh
    if isa(mesh, FEMesh)
        n_total = mesh.number_nodes
        sim.n_cells = mesh.number_cells
    else
        n_total = mesh.nx * mesh.ny * mesh.nz
    end

    sim.n_total = n_total
    sim.spin = create_zeros(3 * n_total)
    sim.prespin = create_zeros(3 * n_total)
    sim.field = create_zeros(3 * n_total)
    sim.energy = create_zeros(n_total)
    sim.pins = create_zeros(Bool, n_total)

    if isa(mesh, FDMesh)
        sim.mu0_Ms = create_zeros(n_total)
    elseif isa(mesh, FEMesh)
        sim.mu0_Ms = zeros(T, sim.n_cells)
        sim.L_mu = create_zeros(3 * n_total)
    else
        sim.mu_s = create_zeros(n_total)
    end

    sim.driver_name = driver
    sim.driver = create_driver(driver, integrator, n_total)
    sim.interactions = []
    sim.save_data = save_data
    saver_name = @sprintf("%s_%s.txt", name, lowercase(driver))
    if driver === "InertialLLG"
        saver_name = @sprintf("%s_illg.txt", name)
    elseif driver === "SpatialLLG"
        saver_name = @sprintf("%s_llg.txt", name)
    end
    sim.saver = create_saver(saver_name, driver)

    if isa(mesh, FDMesh)
        @info "MicroSim (FD) has been created."
    elseif isa(mesh, FEMesh)
        @info "MicroSim (FE) has been created."
    else
        @info "AtomisticSim has been created."
    end

    if global_client != nothing && isa(mesh, FDMesh)
        response = Dict(
            "type" => "fd_mesh_data",
            "success" => success,
            "fd_mesh_data" => Dict(
                "nx" => mesh.nx,
                "ny" => mesh.ny,
                "nz" => mesh.nz,
                "dx" => mesh.dx*1e9,
                "dy" => mesh.dy*1e9,
                "dz" => mesh.dz*1e9,
            )
        )
        send_message(global_client, "run_code_response", response)
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
    set_Ms(sim::AbstractSim, shape::CSGShape, Ms::Number)

Set the saturation magnetization Ms within the Shape.
"""
function set_Ms(sim::AbstractSim, shape::CSGShape, Ms::Number)
    init_scalar!(sim.mu0_Ms, sim.mesh, shape, Ms * mu_0)
    return true
end

"""
    set_pinning(sim::AbstractSim, ids::ArrayOrFunction)

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
function set_pinning(sim::AbstractSim, ids::ArrayOrFunction)
    pins = zeros(Bool, sim.n_total)
    init_scalar!(pins, sim.mesh, ids)
    copyto!(sim.pins, pins)
    return true
end

"""
    set_alpha(sim::AbstractSim, alpha::ArrayOrFunction)

Set the damping parameter α for a simulation with `SpatialLLG` driver.

This function updates the spatially varying damping coefficient array in the simulation's driver.
The α parameter can be specified either as a uniform value (via a function) or as a spatially varying array.

# Arguments
- `sim::AbstractSim`: The simulation object containing the driver
- `alpha::ArrayOrFunction`: The damping parameter specification, which can be:
  - A function `f(i, j, k, dx, dy, dz) -> Float64` that returns the α value at each spatial coordinate
  - An array of Float64 values with length equal to `sim.n_total`

# Examples
```julia
# Set uniform α value of 0.1 using a function
set_alpha(sim, (i, j, k, dx, dy, dz) -> 0.1)

# Set spatially varying α using an array
alpha_array = rand(sim.n_total) * 0.2 + 0.01  # Random values between 0.01 and 0.21
set_alpha(sim, alpha_array)

# Set spatial α with a function
function spatial_alpha(i, j, k, dx, dy, dz)
    if i < 10
        return 0.01
    else
        return 0.99
    end
end
set_alpha(sim, spatial_alpha)
```
"""
function set_alpha(sim::AbstractSim, alpha::ArrayOrFunction)
    if !isa(sim.driver, SpatialLLG)
        @warn("set_alpha only works for the SpatialLLG driver")
        return
    end
    alpha_array = zeros(Float[], sim.n_total)
    init_scalar!(alpha_array, sim.mesh, alpha)
    copyto!(sim.driver.alpha, alpha_array)
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

@doc raw"""
Compute the average magnetization defined as

```math
  \langle \mathbf{m} \rangle  = \frac{1}{V} \int_{V} \mathbf{m} \mathrm{d}V
```
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

    if isdefined(sim.driver, :integrator) && sim.driver.integrator isa AdaptiveRK
        set_initial_condition!(sim, sim.driver.integrator)
    end

    return true
end

"""
    set_driver_arguments(sim::AbstractSim, args::Dict)

Set the parameters of the driver. 
"""
function set_driver_arguments(sim::AbstractSim, args::Dict)

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

    #FIXME: check the type of sim.driver.integrator
    if haskey(args, :tol) && startswith(driver, "LLG")
        sim.driver.integrator.tol = args[:tol]
        delete!(args, :tol)
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
    relax(sim::AbstractSim; max_steps=10000, stopping_dmdt=0.01, save_data_every=0, 
           save_m_every=0, using_time_factor=true, output="ovf", update_time=1.0)

Relaxes the system using either the `LLG` or `SD` driver. This function is compatible with both [Micromagnetic model](@ref) and [Atomistic spin model](@ref).

**Arguments:**

- `max_steps::Int`: Maximum number of steps to run the simulation. Default is `10000`.
- `stopping_dmdt::Float64`: Primary stopping condition for both `LLG` and `SD` drivers. For standard micromagnetic simulations, typical values range from `0.01` to `1`. In `SD` driver mode, where time is not strictly defined, a factor of `γ` is applied to make it comparable to the `LLG` driver. For atomistic models using dimensionless units, set `using_time_factor` to `false` to disable this factor.
- `save_data_every::Int`: Interval for saving overall data such as energies and average magnetization. A negative value disables data saving while `save_data_every=0` saves data only at the end of the relaxation.
- `save_m_every::Int`: Interval for saving magnetization data. A negative value disables magnetization saving while `save_m_every=0` saves magnetization only at the end of the relaxation.
- `using_time_factor::Bool`: Boolean flag to apply a time factor in `SD` mode for comparison with `LLG` mode. Default is `true`.
- `update_time::Float64`: Time interval (in seconds) for updating visualization data. Default is `1.0`.

**Examples:**

```julia
    relax(sim, max_steps=10000, stopping_dmdt=0.1)
```
"""
function relax(sim::AbstractSim; max_steps=10000, stopping_dmdt=0.01, save_data_every=0,
               save_m_every=0, using_time_factor=true, maxsteps=nothing, output="ovf", update_time=1.0)
    if !isnothing(maxsteps)
        @warn "The parameter 'maxsteps' is deprecated. Please use 'max_steps' in relax function."
        max_steps = maxsteps
    end

    # to dertermine which driver is used.
    llg_driver = isa(sim.driver, LLG)

    time_factor = using_time_factor ? 2.21e5 / 2 : 1.0
    dmdt_factor = using_time_factor ? (2 * pi / 360) * 1e9 : 1

    output_folder = @sprintf("%s_%s", sim.name, nameof(typeof(sim.driver)))
    if save_m_every>0
        mkpath(output_folder)
    end

    N_spins = sim.n_total
    dm = create_zeros(3 * N_spins)
    
    last_update_time = time()

    driver = sim.driver
    @info @sprintf("Running Driver : %s.", typeof(driver))
    for i in 1:max_steps
        @timeit timer "run_step" run_step(sim, driver)

        if llg_driver && isa(sim.driver.integrator, AdaptiveRK)
            copyto!(sim.prespin, sim.spin)
            y = view(sim.driver.integrator.y_current, 1:(3 * sim.n_total))
            copyto!(sim.spin, y)
        end

        step_size = llg_driver ? driver.integrator.step : driver.tau / time_factor

        compute_dm!(dm, sim.prespin, sim.spin, N_spins)
        max_dmdt = maximum(dm) / step_size

        stop_flag = max_dmdt < stopping_dmdt * dmdt_factor

        t = llg_driver ? sim.driver.integrator.t : 0.0
        if llg_driver
            Verbose[] &&
                @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
                               i, step_size, t, max_dmdt / dmdt_factor)
        else
            Verbose[] &&
                @info @sprintf("step =%5d  step_size=%10.6e  max_dmdt=%10.6e", i, step_size,
                               max_dmdt / dmdt_factor)
        end

        if sim.save_data
            sim.saver.nsteps += 1
        end

        if (save_data_every > 0 && i % save_data_every == 0) || (stop_flag && save_data_every == 0)
            compute_system_energy(sim, sim.spin, t)
            sim.saver.t = t
            write_data(sim)
        end

        if (save_m_every > 0 && i % save_m_every == 0) || (stop_flag && save_m_every == 0)
            if output == "ovf"
                save_ovf(sim, joinpath(output_folder, @sprintf("m_%08d.ovf", i)))
            elseif output == "vts" || output == "vtu"
                save_vtk_points(sim, joinpath(output_folder, @sprintf("m_%08d", i)))
            end
        end

        current_time = time()
        if update_time > 0 && (current_time - last_update_time >= update_time || stop_flag)
            if global_client != nothing
                response = Dict(
                    "type" => "m_data",
                    "m_data" => Array(sim.spin),
                )
                send_message(global_client, "m_data", response)
                last_update_time = current_time
            end
        end

        if stop_flag
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g @steps=%g, Done!",
                           stopping_dmdt, i)
            break
        end

    end

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
    if save_data
        compute_system_energy(sim, sim.spin, t_end)
        write_data(sim)
    end
    sim.saver.nsteps += 1

    return nothing
end

function run_until(sim::AbstractSim, t_end::Float64, integrator::Integrator,
                   save_data::Bool)
    if t_end < integrator.t - integrator.step
        @info("Run_until: t_end >= integrator.t - integrator.step")
        return
    elseif t_end == integrator.t
        sim.saver.t = t_end
        if save_data
            compute_system_energy(sim, sim.spin, t_end)
            write_data(sim)
        end
        sim.saver.nsteps += 1
        return
    end

    copyto!(sim.prespin, sim.spin)
    integrate_to_time(sim, integrator, t_end)
    y = view(integrator.y_current, 1:(3 * sim.n_total))
    copyto!(sim.spin, y)

    sim.saver.t = t_end
    if save_data
        compute_system_energy(sim, sim.spin, t_end)
        write_data(sim)
    end
    sim.saver.nsteps += 1
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
        sim.driver = create_driver(driver, integrator, sim.n_total)
        sim.driver_name = driver
        if sim.driver.integrator isa AdaptiveRK
            set_initial_condition!(sim, sim.driver.integrator)
        end

        saver = sim.saver

        saver_name = @sprintf("%s_%s.txt", sim.name, lowercase(driver))
        if driver === "InertialLLG"
            saver_name = @sprintf("%s_illg.txt", sim.name)
        elseif driver === "SpatialLLG"
            saver_name = @sprintf("%s_llg.txt", sim.name)
        end

        saver.name = saver_name
        saver.t = 0
        saver.nsteps = 0
        saver.header_saved = false
        contains_time = false
        for item in saver.items
            if item.name == "time"
                contains_time = true
            end
        end
        if driver != "SD" && !contains_time
            time = SaverItem("time", "<s>", o::AbstractSim -> o.saver.t)
            insert!(saver.items, 2, time)
            #push!(saver.items, time)
        end
    end

    return set_driver_arguments(sim, args)
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
- `dmi_type` : the type of DMI, could be "bulk", "interfacial" or "D2d".
- `Ku`: the anisotropy constant, should be [`NumberOrArrayOrFunction`](@ref).
- `axis`: the anisotropy axis, should be a tuple, such as (0,0, 1)
- `Kc`: the cubic anisotropy constant, should be [`NumberOrArrayOrFunction`](@ref).
- `axis1`: the cubic anisotropy axis1, should be a tuple, such as (1,0,0)
- `axis2`: the cubic anisotropy axis2, should be a tuple, such as (0,1,0)
- `demag` : include demagnetization or not, should be a boolean, i.e., true or false. By default,  demag=false.
- `H`: the external field, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref).
- `m0` : the initial magnetization, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref).
- `T` : the temperature, should be should be [`NumberOrArrayOrFunction`](@ref).
- `shape` : the shape defines the geometry of the sample, where parameters are configured.
"""
function create_sim(mesh; args...)
    return create_sim(mesh, Dict(args))
end

function create_sim(mesh, args::Dict)
    #Ms=8e5, A=1.3e-11, D=0, Ku=0, axis=(0,0,1), H=(0,0,0), m0=(0,0,1), demag=false, driver="SD", name="relax"
    driver = haskey(args, :driver) ? args[:driver] : "SD"
    name = haskey(args, :name) ? args[:name] : "unnamed"
    shape = haskey(args, :shape) ? args[:shape] : nothing

    integrator = get(args, :integrator, "DormandPrince")

    #Create the mesh using given driver and name
    sim = Sim(mesh; driver=driver, integrator=integrator, name=name)

    #If the simulation is the standard micromagnetics simulation.
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

        for key in [:Ms, :A, :D, :demag, :dmi_type]
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

    if haskey(args, :Kc)
        axis1 = haskey(args, :axis1) ? args[:axis1] : (1, 0, 0)
        axis2 = haskey(args, :axis2) ? args[:axis2] : (0, 1, 0)

        add_cubic_anis(sim, args[:Kc]; axis1=axis1, axis2=axis2)

        haskey(args, :axis1) && delete!(args, :axis1)
        haskey(args, :axis2) && delete!(args, :axis2)
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
    set_driver_arguments(sim, args)

    # set m0 anyway
    m0_value = haskey(args, :m0) ? args[:m0] : (0.8, 0.6, 0)
    init_m0(sim, m0_value)

    for key in [:driver, :name, :Ku, :Kc, :H, :m0, :shape, :T]
        haskey(args, key) && delete!(args, key)
    end

    for key in args
        #@warn @sprintf("Key '%s' is not used.", key)
    end

    return sim
end

"""
    run_sim(sim::AbstractSim; steps=10, dt=1e-12, save_data=true, save_m_every=1, saver_item=nothing, call_back=nothing)

Run the simulation for a total duration of `steps * dt` seconds.

- `steps`: The total number of simulation steps. The total simulation time will be `steps * dt`.
- `dt`: The time step for each simulation step (in seconds). It controls the time interval between successive steps.
- `save_data`: A boolean flag to control whether the overall simulation data (such as total energy, exchange energy, and average magnetization) is saved at each step. The default is `true`.
- `save_m_every`: Specifies how often to save the magnetization configuration. For example, if `save_m_every = 1`, the magnetization is saved at every step. A negative value will disable magnetization saving entirely.
- `saver_item`: A `SaverItem` instance or a list of `SaverItem` instances. These are custom data-saving utilities that can be used to store additional quantities during the simulation (e.g., guiding centers or other derived values). If `nothing`, no additional data is saved beyond the default.
- `call_back`: A user-defined function or `nothing`. If provided, this function will be called at every step, allowing for real-time inspection or manipulation of the simulation state.
- `update_time`: The time interval (in seconds) at which the magnetization data is sent to the client. The default is `1.0` second.

#### Example Usage

To compute a custom quantity, such as the guiding center of the magnetization, and save it in the output table, you can define a `SaverItem` as shown below:

```julia
run_sim(sim, saver_item=SaverItem(("Rx", "Ry"), ("<m>", "<m>"), compute_guiding_center))
```

In this example:
- `("Rx", "Ry")` are the labels for the data in the output.
- `("<m>", "<m>")` are the units of guiding center.
- `compute_guiding_center` is the function that computes the guiding center at each step.

This setup will save the computed guiding center to the simulation output, in addition to the default data like energies and average magnetization.
"""
function run_sim(sim::AbstractSim; steps=10, dt=1e-10, save_data=true, save_m_every=-1,
                 saver_item=nothing, call_back=nothing, output="ovf", update_time=1.0)
    if isa(saver_item, SaverItem)
        push!(sim.saver.items, saver_item)
    end

    if isa(saver_item, AbstractArray)
        for item in saver_item
            push!(sim.saver.items, item)
        end
    end

    output_folder = @sprintf("%s_%s", sim.name, nameof(typeof(sim.driver)))
    if save_m_every >= 0
        mkpath(output_folder)
    end

    last_update_time = time()
    for i in 0:steps
        run_until(sim, i * dt; save_data=save_data)

        !Verbose[] && (steps == i ? println() : print("."))
        (Verbose[] || steps == i) && @info @sprintf("step =%5d  t = %10.6e", i, i * dt)

        call_back !== nothing && call_back(sim, i * dt)

        if save_m_every > 0 && (i - 1) % save_m_every == 0
            if output == "ovf"
                save_ovf(sim, joinpath(output_folder, @sprintf("m_%08d.ovf", i)))
            elseif output == "vts" || output == "vtu"
                save_vtk_points(sim, joinpath(output_folder, @sprintf("m_%08d", i)))
            end
        end

        current_time = time()
        if update_time > 0 && (current_time - last_update_time >= update_time || i == steps)
            if global_client != nothing
                response = Dict(
                    "type" => "m_data",
                    "m_data" => Array(sim.spin),
                )
                send_message(global_client, "m_data", response)
                last_update_time = current_time
            end
        end
    end

    return nothing
end

"""
    sim_with(args::Union{NamedTuple, Dict})

[High-Level Interface](@ref) for starting a typical micromagnetic simulation. All parameters are set using `args`, which can be either a `NamedTuple` or a `Dict`.

# Keywords

- `mesh`: A mesh must be provided to start the simulation. The mesh could be [`FDMesh`](@ref), [`CubicMesh`](@ref), or [`TriangularMesh`](@ref).
- `name`: The name of the simulation, provided as a string.
- `task`: The type of simulation task, which can be `"Relax"` or `"Dynamics"`. The default is "Relax".
- `driver`: The name of the driver, which should be "SD", "LLG", or "LLG_STT". The default is "SD".
- `alpha`: The Gilbert damping parameter in the LLG equation, provided as a number.
- `beta`: The nonadiabatic strength in the LLG equation with spin-transfer torques (Zhang-Li model), provided as a number.
- `gamma`: The gyromagnetic ratio, with a default value of 2.21e5.
- `ux`, `uy`, `uz`: The components of the spin-transfer torque strength.
- `ufun`: A time-dependent function for `u`.
- `Ms`: The saturation magnetization, which should be a [`NumberOrArrayOrFunction`](@ref). The default is `Ms=8e5`.
- `mu_s`: The magnetic moment, which should be a [`NumberOrArrayOrFunction`](@ref). The default is `mu_s=2*mu_B`.
- `A` or `J`: The exchange constant, which should be a [`NumberOrArrayOrFunction`](@ref).
- `D`: The Dzyaloshinskii-Moriya interaction (DMI) constant, which should be a [`NumberOrArrayOrFunction`](@ref).
- `dmi_type`: The type of DMI, either "bulk" or "interfacial".
- `Ku`: The anisotropy constant, which should be a [`NumberOrArrayOrFunction`](@ref).
- `axis`: The anisotropy axis, provided as a tuple, e.g., `(0, 0, 1)`.
- `Kc`: the cubic anisotropy constant, should be [`NumberOrArrayOrFunction`](@ref).
- `axis1`: the cubic anisotropy axis1, should be a tuple, such as (1,0,0)
- `axis2`: the cubic anisotropy axis2, should be a tuple, such as (0,1,0)
- `demag`: Whether to include demagnetization. This should be a boolean (`true` or `false`). The default is `demag=false`.
- `H`: The external magnetic field, which should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref).
- `m0`: The initial magnetization, which should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref).
- `T`: The temperature, which should be a [`NumberOrArrayOrFunction`](@ref).
- `shape`: The shape defines the geometry of the sample, where parameters are configured.
- `steps`: The total number of simulation steps for the `Dynamics` task.
- `dt`: The time interval of each step, so the total simulation time is `steps * dt` for the `Dynamics` task.
- `max_steps::Int`: Maximum number of steps to run the simulation for the `Relax` task. Default is `10000`.
- `saver_item`: A `SaverItem` instance or a list of `SaverItem` instances. These are custom data-saving utilities that can be used to store additional quantities during the simulation (e.g., guiding centers or other derived values). If `nothing`, no additional data is saved beyond the default.
- `call_back`: A user-defined function or `nothing`. If provided, this function will be called at every step, allowing for real-time inspection or manipulation of the simulation state.
- `stopping_dmdt::Float64`: Primary stopping condition for both `LLG` and `SD` drivers. For standard micromagnetic simulations, typical values range from `0.01` to `1`. In `SD` driver mode, where time is not strictly defined, a factor of `γ` is applied to make it comparable to the `LLG` driver. For atomistic models using dimensionless units, set `using_time_factor` to `false` to disable this factor.
- `relax_data_every::Int`: Interval for saving overall data such as energies and average magnetization during a `Relax` task. A negative value disables data saving (e.g., `relax_data_every = -1` saves data only at the end of the relaxation).
- `dynamic_data_save::Bool`: Boolean flag to enable or disable saving overall data such as energies and average magnetization during the `Dynamics` task. Set to `true` to enable, or `false` to disable.
- `relax_m_every::Int`: Interval for saving magnetization data during a `Relax` task. A negative value disables magnetization saving while `relax_m_every=0` saves magnetization only at the end of the relaxation.
- `dynamic_m_every ::Int`: Interval for saving magnetization data during a `Dynamics` task. A negative value disables magnetization saving while `dynamic_m_every =0` saves magnetization only at the end of the dynamics.
- `using_time_factor::Bool`: Boolean flag to apply a time factor in `SD` mode for comparison with `LLG` mode. Default is `true`.
- `save_vtk::Bool`: Boolean flag to save the magnetization to vtk files after finishing each task. Default is `false`.

#### Example

See examples at [High-Level Interface](@ref).

#### Notes

- **Suffix Usage**: Parameters like `H`, `Ms`, `Ku`, `A`, `D`, `task`, and `driver` can be defined as arrays using the `_s` or `_sweep` suffix. When using these suffixes, ensure that the corresponding array lengths match. This allows the simulation to iterate over different values for these parameters.
- **Argument Types**: The `args` parameter can be either a `NamedTuple` or a `Dict`, providing flexibility in how you organize and pass the simulation parameters.
- **Driver Selection**: The `driver` parameter (or `driver_s` for multiple drivers) specifies the simulation type. Options include `"SD"` for the steepest-descent method, `"LLG"` for the Landau-Lifshitz-Gilbert equation, and `"LLG_STT"` for simulations involving spin-transfer torques.
- **Stopping Criterion**: The `stopping_dmdt` parameter is critical for determining when to stop a simulation, particularly in relaxation tasks. It measures the rate of change in magnetization, with typical values ranging from `0.01` to `1`. For atomistic models, the `using_time_factor` flag can be set to `false` to disable time scaling.
- **Data Saving**: The `relax_m_every` and `dynamic_m_every ` parameters control how frequently magnetization data is saved during `Relax` and `Dynamics` tasks, respectively. Use negative values to disable data saving for these tasks.

"""
function sim_with(args::Union{NamedTuple,Dict})
    #convert args to a dict
    if !isa(args, Dict)
        args = Dict(key => value for (key, value) in pairs(args))
    end

    task = get(args, :task, "Relax")
    driver = get(args, :driver, "SD")
    stopping_dmdt = get(args, :stopping_dmdt, 0.1)
    max_steps = get(args, :max_steps, 10000)
    relax_data_every = get(args, :relax_data_every, 0)
    relax_m_every = get(args, :relax_m_every, 0)
    using_time_factor = get(args, :using_time_factor, true)
    steps = get(args, :steps, 100)
    dt = get(args, :dt, 1e-11)
    dynamic_data_save = get(args, :dynamic_data_save, true)
    dynamic_m_every  = get(args, :dynamic_m_every , -1)
    call_back = get(args, :call_back, nothing)
    saver_item = get(args, :saver_item, nothing)
    vtk_saving = get(args, :save_vtk, false)
    name = haskey(args, :name) ? args[:name] : "unnamed"

    if !haskey(args, :mesh)
        error("A mesh must be provided to start the simulation.")
    end

    mesh = args[:mesh]
    delete!(args, :mesh)

    N = check_sweep_lengths(args)

    # common single task without range
    if N == 0
        sim = create_sim(mesh, args)

        haskey(args, :task) && delete!(args, :task)

        if isa(saver_item, SaverItem)
            push!(sim.saver.items, saver_item)
        end

        if isa(saver_item, AbstractArray)
            for item in saver_item
                push!(sim.saver.items, item)
            end
        end

        if startswith(lowercase(task), "rel") # task == relax 
            for key in [:stopping_dmdt, :max_steps, :relax_data_every, :relax_m_every,
                        :using_time_factor]
                haskey(args, key) && delete!(args, key)
            end

            relax(sim; max_steps=max_steps, stopping_dmdt=stopping_dmdt,
                  save_data_every=relax_data_every, save_m_every=relax_m_every,
                  using_time_factor=using_time_factor)

        elseif startswith(lowercase(task), "dyn") # task == "dynamics"
            if driver == "SD"
                set_driver(sim; driver="LLG", integrator="DormandPrince")
                set_driver_arguments(sim, args)
                set_initial_condition!(sim, sim.driver.integrator)
            end

            for key in [:steps, :dt, :dynamic_data_save, :dynamic_m_every , :call_back,
                        :saver_item, :save_vtk]
                haskey(args, key) && delete!(args, key)
            end

            run_sim(sim; steps=steps, dt=dt, save_data=dynamic_data_save,
                    save_m_every=dynamic_m_every , call_back=call_back)
        end

        for key in args
            @warn @sprintf("Key '%s' is not used.", key)
        end

        if vtk_saving
            !isdir("vtks") && mkdir("vtks")
            vtkname = @sprintf("vtks/%s.vts", name)
            save_vtk(sim, vtkname)
        end

        return sim
    end

    # Now we need to deal with the case that N > 0
    dict = extract_sweep_keys(args)
    for (k, v) in dict
        args[k] = v[1]
    end

    shape = get(args, :shape, nothing)

    sim = create_sim(mesh, args)
    haskey(args, :task) && delete!(args, :task)

    if isa(saver_item, SaverItem)
        push!(sim.saver.items, saver_item)
    end

    if isa(saver_item, AbstractArray)
        for item in saver_item
            push!(sim.saver.items, item)
        end
    end

    for n in 1:N
        task_ = haskey(dict, :task) ? dict[:task][n] : task
        driver_ = haskey(dict, :driver) ? dict[:driver][n] : driver

        if haskey(dict, :Ms)
            Ms_ = dict[:Ms][n]
            if shape === nothing
                set_Ms(sim, Ms_)
            else
                set_Ms(sim, shape, Ms_)
            end
        end

        if haskey(dict, :H)
            update_zeeman(sim, dict[:H][n])
        end

        # TODO: we can add more ...

        if startswith(lowercase(task_), "rel")
            for key in [:stopping_dmdt, :max_steps, :relax_data_every, :relax_m_every,
                        :using_time_factor]
                haskey(args, key) && delete!(args, key)
            end
            relax(sim; max_steps=max_steps, stopping_dmdt=stopping_dmdt,
                  save_data_every=relax_data_every, save_m_every=relax_m_every,
                  using_time_factor=using_time_factor)

        elseif startswith(lowercase(task_), "dyn")
            if driver_ == "SD"
                set_driver(sim; driver="LLG", integrator="DormandPrince")
            else
                set_driver(sim; driver=driver_, integrator="DormandPrince")
            end
            set_driver_arguments(sim, args)
            set_initial_condition!(sim, sim.driver.integrator)

            for key in [:steps, :dt, :dynamic_data_save, :dynamic_m_every , :call_back,
                        :saver_item, :save_vtk]
                haskey(args, key) && delete!(args, key)
            end

            run_sim(sim; steps=steps, dt=dt, save_data=dynamic_data_save,
                    save_m_every=dynamic_m_every , call_back=call_back)

        else
            error("Only support two types of task: 'Relax' and 'Dynamics'.")
        end

        if vtk_saving
            !isdir("vtks") && mkdir("vtks")
            vtkname = @sprintf("vtks/%s_%d.vts", name, n)
            save_vtk(sim, vtkname)
        end
    end

    for key in args
        @warn @sprintf("Key '%s' is not used.", key)
    end

    return sim
end

function sim_with(; args...)
    return sim_with(Dict(args))
end

function advance_step(sim::AbstractSim)
    return advance_step(sim, sim.driver.integrator)
end









