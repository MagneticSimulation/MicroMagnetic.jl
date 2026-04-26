export create_sim, sim_with

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