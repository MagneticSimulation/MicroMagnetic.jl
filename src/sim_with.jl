export sim_with

# --------------------------------------------------------------------------
# Keyword vocabularies used for the fail-fast validation of `sim_with` (I-09).
# Every keyword accepted by `sim_with` must be listed here; unknown keys throw
# an error before any simulation is started.
# --------------------------------------------------------------------------
const _COMMON_KEYS = (:mesh, :name, :task, :driver, :integrator, :save_vtk, :saver_item,
                      :stt, :sot)
const _MATERIAL_KEYS = (:Ms, :mu_s, :A, :J, :D, :dmi_type, :Ku, :axis, :Kc, :axis1,
                        :axis2, :demag, :H, :T, :m0, :shape)
const _DRIVER_KEYS = (:alpha, :gamma, :tol)
# keys removed together with the LLG_STT/LLG_CPP drivers in v0.6.0 (I-10)
const _REMOVED_KEYS = Dict(:beta => "put the nonadiabatic parameter inside `stt=(model=:zhang_li, ..., xi=...)`",
                           :ux => "use `stt=(model=:zhang_li, b=<u>, J=<direction>, ...)`",
                           :uy => "use `stt=(model=:zhang_li, b=<u>, J=<direction>, ...)`",
                           :uz => "use `stt=(model=:zhang_li, b=<u>, J=<direction>, ...)`",
                           :ufun => "use the time dependence inside `stt=(..., ft=<fun>)`")
const _RELAX_KEYS = (:stopping_dmdt, :max_steps, :relax_data_every, :relax_m_every,
                     :using_time_factor)
const _DYN_KEYS = (:steps, :dt, :dynamic_data_save, :dynamic_m_every, :call_back)
const _SIM_WITH_KEYS = (_COMMON_KEYS..., _MATERIAL_KEYS..., _DRIVER_KEYS...,
                        _RELAX_KEYS..., _DYN_KEYS...)

# The only parameters that support per-stage sweeping via `_s`/`_sweep` (I-07).
const _SWEEP_KEYS = (:task, :driver, :Ms, :H)

_is_relax(task) = startswith(lowercase(string(task)), "rel")
_is_dynamics(task) = startswith(lowercase(string(task)), "dyn")
_is_llg_driver(name) = startswith(lowercase(string(name)), "llg")

function _check_task(task)
    if !_is_relax(task) && !_is_dynamics(task)
        error("`sim_with` only supports task = \"Relax\" or \"Dynamics\", got `$task`.")
    end
    return nothing
end

function _levenshtein(a::AbstractString, b::AbstractString)
    n, m = length(a), length(b)
    (n == 0 || m == 0) && return max(n, m)
    prev = collect(0:m)
    cur = similar(prev)
    for i in 1:n
        cur[1] = i
        ca = a[i]
        for j in 1:m
            cost = ca == b[j] ? 0 : 1
            cur[j + 1] = min(cur[j] + 1, prev[j + 1] + 1, prev[j] + cost)
        end
        prev, cur = cur, prev
    end
    return prev[m + 1]
end

function _suggest_key(key::Symbol)
    best, bestd = :task, typemax(Int)
    for cand in _SIM_WITH_KEYS
        d = _levenshtein(String(key), String(cand))
        if d < bestd
            best, bestd = cand, d
        end
    end
    return bestd <= 3 ? " (did you mean `:$best`?)" : ""
end

# A key is known when it is listed in `_SIM_WITH_KEYS` or when it is the `_s`/`_sweep`
# variant of a known key (the per-stage support of the sweep variant is checked later, I-07).
function _is_known_key(key::Symbol)
    key in _SIM_WITH_KEYS && return true
    name = String(key)
    if endswith(name, "_sweep")
        return Symbol(chop(name; tail = 6)) in _SIM_WITH_KEYS
    elseif endswith(name, "_s") && name != "mu_s"
        return Symbol(chop(name; tail = 2)) in _SIM_WITH_KEYS
    end
    return false
end

function _check_known_keys(args::Dict)
    for key in keys(args)
        if haskey(_REMOVED_KEYS, key)
            error("`sim_with` no longer accepts `:$key` (removed in v0.6.0): " *
                  _REMOVED_KEYS[key] * ".")
        end
        if !_is_known_key(key)
            error("`sim_with` got an unknown key `:$key`$(_suggest_key(key)). " *
                  "See `?sim_with` for the list of supported keys.")
        end
    end
    return nothing
end

# Warn (before running) about keys that cannot apply to the tasks being run.
function _warn_unused_keys(args::Dict, plausible::Tuple)
    for key in sort!(collect(keys(args)); by = string)
        key in plausible || @warn "Key `:$key` is not used."
    end
    return nothing
end

# Assemble a Sim from a mesh plus a material/initial-state keyword Dict. This is
# a private helper of `sim_with` (the public entry point); the keyword reference
# lives in the `sim_with` docstring.
function _build_sim(mesh, args::Dict)
    # never modify the caller's Dict (I-08)
    args = copy(args)

    # warn about orphaned keys that will not be used (I-14)
    haskey(args, :dmi_type) && !haskey(args, :D) &&
        @warn "Key `:dmi_type` is not used because `:D` is not provided."
    haskey(args, :axis) && !haskey(args, :Ku) &&
        @warn "Key `:axis` is not used because `:Ku` is not provided."
    haskey(args, :axis1) && !haskey(args, :Kc) &&
        @warn "Key `:axis1` is not used because `:Kc` is not provided."
    haskey(args, :axis2) && !haskey(args, :Kc) &&
        @warn "Key `:axis2` is not used because `:Kc` is not provided."

    driver = haskey(args, :driver) ? args[:driver] : "SD"
    name = haskey(args, :name) ? args[:name] : "unnamed"
    shape = haskey(args, :shape) ? args[:shape] : nothing

    integrator = get(args, :integrator, "DormandPrince")

    #Create the mesh using given driver and name
    sim = Sim(mesh; driver=driver, integrator=integrator, name=name)

    # Spin-transfer / spin-orbit torques are added as interactions; they only make
    # sense for LLG-family drivers because the non-conservative torque field must
    # not enter the SD energy minimization (I-10)
    if haskey(args, :stt)
        if _is_llg_driver(driver)
            add_stt(sim; args[:stt]...)
        else
            @warn "Key `:stt` is ignored because driver=\"$driver\" is not an LLG-family driver."
        end
        delete!(args, :stt)
    end
    if haskey(args, :sot)
        if _is_llg_driver(driver)
            add_sot(sim; args[:sot]...)
        else
            @warn "Key `:sot` is ignored because driver=\"$driver\" is not an LLG-family driver."
        end
        delete!(args, :sot)
    end

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

        # warn about material keys ignored by this mesh type (I-14)
        haskey(args, :mu_s) &&
            @warn "Key `:mu_s` is ignored for FDMesh (use `:Ms` instead)."
        haskey(args, :J) &&
            @warn "Key `:J` is ignored for FDMesh (use `:A` instead)."

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

        # warn about material keys ignored by this mesh type (I-14)
        haskey(args, :Ms) &&
            @warn "Key `:Ms` is ignored for atomistic meshes (use `:mu_s` instead)."
        haskey(args, :A) &&
            @warn "Key `:A` is ignored for atomistic meshes (use `:J` instead)."
    else
        error("Unsupported mesh type $(typeof(mesh)). Only FDMesh and AtomisticMesh " *
              "(CubicMesh, TriangularMesh, ...) are supported by `_build_sim`.")
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

    for key in [:driver, :name, :integrator, :Ku, :Kc, :H, :m0, :shape, :T]
        haskey(args, key) && delete!(args, key)
    end

    return sim
end

"""
    sim_with(args::Union{NamedTuple, Dict})

[High-Level Interface](@ref) for starting a typical micromagnetic simulation. All parameters are set using `args`, which can be either a `NamedTuple` or a `Dict`.

# Keywords

- `mesh`: A mesh must be provided to start the simulation. The mesh could be [`FDMesh`](@ref), [`CubicMesh`](@ref), or [`TriangularMesh`](@ref).
- `name`: The name of the simulation, provided as a string. Default is "unnamed".
- `task`: The type of simulation task, which can be `"Relax"` or `"Dynamics"`. The default is "Relax".
- `driver`: The name of the driver, which should be "SD" or "LLG". The default is "SD".
- `integrator`: The integrator used to create the simulation, e.g. `"DormandPrince"`. Default is `"DormandPrince"`.
- `alpha`: The Gilbert damping parameter in the LLG equation, provided as a number.
- `gamma`: The gyromagnetic ratio, with a default value of 2.21e5.
- `stt`: Parameters forwarded to [`add_stt`](@ref) as a `NamedTuple`, e.g. `stt=(model=:zhang_li, b=-72.438, J=(1,0,0), xi=0.05)`. The torque is applied only in stages that run with an LLG-family driver; SD stages are not affected.
- `sot`: Parameters forwarded to [`add_sot`](@ref) as a `NamedTuple`, applied under the same driver rule as `stt`.
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
- `relax_data_every::Int`: Interval for saving overall data such as energies and average magnetization during a `Relax` task. `0` saves data only at the end of the relaxation and a negative value disables data saving.
- `relax_m_every::Int`: Interval for saving magnetization data during a `Relax` task. `0` saves magnetization only at the end of the relaxation (also when `max_steps` is reached without convergence) and a negative value disables magnetization saving.
- `dynamic_data_save::Bool`: Boolean flag to enable or disable saving overall data such as energies and average magnetization during the `Dynamics` task. Set to `true` to enable, or `false` to disable.
- `dynamic_m_every::Int`: Interval for saving magnetization data during a `Dynamics` task. `0` saves magnetization only at the end of the dynamics and a negative value disables magnetization saving.
- `using_time_factor::Bool`: Boolean flag to apply a time factor in `SD` mode for comparison with `LLG` mode. Default is `true`.
- `save_vtk::Bool`: Boolean flag to save the magnetization to vtk files after finishing each task. Default is `false`.

#### Example

See examples at [High-Level Interface](@ref).

#### Notes

- **Fail-fast validation**: Unknown keys throw an error immediately (before the simulation starts), with a "did you mean" hint for likely typos. Keys that cannot apply to the configured task(s) produce a warning before the simulation runs.
- **Argument Types**: The `args` parameter can be either a `NamedTuple` or a `Dict`. Neither is modified; `sim_with` works on an internal copy.
- **Suffix Usage**: Only `task`, `driver`, `Ms` and `H` support the `_s` or `_sweep` suffix to iterate over per-stage values (e.g. `task_s`, `H_s`). The corresponding array lengths must match. Sweeping over any other key throws an error.
- **Driver Selection**: The `driver` parameter (or `driver_s` for multiple stages) specifies the simulation type. Options include `"SD"` for the steepest-descent method and `"LLG"` for the Landau-Lifshitz-Gilbert equation. The `"LLG_STT"`/`"LLG_CPP"` drivers were removed in v0.6.0; use `stt=(...)`/`sot=(...)` together with the `LLG` driver instead. `driver_s` is honored for both `Relax` and `Dynamics` stages, and driver parameters (`alpha`, `gamma`, ...) are re-applied whenever the driver changes.
- **Stopping Criterion**: The `stopping_dmdt` parameter is critical for determining when to stop a simulation, particularly in relaxation tasks. It measures the rate of change in magnetization, with typical values ranging from `0.01` to `1`. For atomistic models, the `using_time_factor` flag can be set to `false` to disable time scaling.
- **Data Saving**: `relax_data_every`/`relax_m_every` and `dynamic_m_every` control how often data/magnetization are saved during `Relax` and `Dynamics` tasks. A value of `0` saves only at the end of the task, a negative value disables saving.

"""
function sim_with(args::Union{NamedTuple,Dict})
    # work on a copy so the caller's Dict/NamedTuple is never modified (I-08)
    args = isa(args, Dict) ? copy(args) :
           Dict(key => value for (key, value) in pairs(args))

    # fail fast on unknown keys, before anything is built or run (I-09)
    _check_known_keys(args)

    haskey(args, :mesh) || error("A mesh must be provided to start the simulation.")
    mesh = pop!(args, :mesh)

    # validate the task value up front (I-04)
    task = get(args, :task, "Relax")
    _check_task(task)
    delete!(args, :task)

    driver = get(args, :driver, "SD")
    integ = get(args, :integrator, "DormandPrince")
    stopping_dmdt = get(args, :stopping_dmdt, 0.1)
    max_steps = get(args, :max_steps, 10000)
    relax_data_every = get(args, :relax_data_every, 0)
    relax_m_every = get(args, :relax_m_every, 0)
    using_time_factor = get(args, :using_time_factor, true)
    steps = get(args, :steps, 100)
    dt = get(args, :dt, 1e-11)
    dynamic_data_save = get(args, :dynamic_data_save, true)
    dynamic_m_every = get(args, :dynamic_m_every, -1)
    call_back = get(args, :call_back, nothing)
    saver_item = get(args, :saver_item, nothing)
    stt_params = pop!(args, :stt, nothing)
    sot_params = pop!(args, :sot, nothing)
    vtk_saving = get(args, :save_vtk, false)
    name = get(args, :name, "unnamed")

    N = check_sweep_lengths(args)
    dict = Dict{Symbol,Any}()
    if N > 0
        dict = extract_sweep_keys(args)
        # only task/driver/Ms/H support per-stage sweeping (I-07)
        for key in keys(dict)
            if !(key in _SWEEP_KEYS)
                error("`sim_with` does not support sweeping over `:$key`. " *
                      "Supported per-stage keys are: :task, :driver, :Ms, :H.")
            end
        end
        for (key, value) in dict
            args[key] = value[1]
        end
    end

    # snapshot of the driver parameters; re-applied to every driver created or
    # switched to below, so switching drivers never falls back to the defaults (I-03)
    driver_kw = Dict(key => args[key] for key in _DRIVER_KEYS if haskey(args, key))

    shape = get(args, :shape, nothing)
    sim = _build_sim(mesh, args)

    # _build_sim works on its own copy, so drop the material keys here; it has
    # either consumed them or warned about the ones it ignores (I-08/I-14)
    for key in _MATERIAL_KEYS
        delete!(args, key)
    end

    # attach custom saver items (both tasks)
    if isa(saver_item, SaverItem)
        push!(sim.saver.items, saver_item)
    elseif isa(saver_item, AbstractArray)
        for item in saver_item
            push!(sim.saver.items, item)
        end
    end
    delete!(args, :saver_item)

    # STT/SOT are applied as soon as a stage runs with an LLG-family driver; SD
    # stages never see the non-conservative torque field (I-10)
    stt_applied = Ref(false)
    sot_applied = Ref(false)
    function _apply_torques!()
        if _is_llg_driver(sim.driver_name)
            if stt_params !== nothing && !stt_applied[]
                add_stt(sim; stt_params...)
                stt_applied[] = true
            end
            if sot_params !== nothing && !sot_applied[]
                add_sot(sim; sot_params...)
                sot_applied[] = true
            end
        end
        return nothing
    end
    _apply_torques!()

    # warn early about keys that cannot apply to the configured tasks (I-09)
    n_stages = max(N, 1)
    stage_task(i) = N > 0 && haskey(dict, :task) ? dict[:task][i] : task
    stage_driver(i) = N > 0 && haskey(dict, :driver) ? dict[:driver][i] : driver
    has_relax = any(i -> _is_relax(stage_task(i)), 1:n_stages)
    has_dyn = any(i -> _is_dynamics(stage_task(i)), 1:n_stages)
    plausible = (_COMMON_KEYS..., _DRIVER_KEYS...,
                 (has_relax ? _RELAX_KEYS : ())...,
                 (has_dyn ? _DYN_KEYS : ())...)
    _warn_unused_keys(args, plausible)
    if (stt_params !== nothing || sot_params !== nothing) &&
       !any(i -> _is_dynamics(stage_task(i)) || _is_llg_driver(stage_driver(i)), 1:n_stages)
        @warn "Keys `:stt`/`:sot` are provided but will not be applied: no stage " *
              "runs with an LLG-family driver."
    end

    # common single task without range
    if N == 0
        if _is_relax(task)
            relax(sim; max_steps=max_steps, stopping_dmdt=stopping_dmdt,
                  save_data_every=relax_data_every, save_m_every=relax_m_every,
                  using_time_factor=using_time_factor)
        else  # dynamics
            if lowercase(string(driver)) == "sd"
                # SD cannot run dynamics; fall back to LLG and re-apply the driver
                # parameters to the freshly created driver (I-03)
                @warn "driver=\"SD\" cannot run a Dynamics stage; switched to \"LLG\"."
                set_driver(sim; driver="LLG", integrator=integ, driver_kw...)
            end
            set_initial_condition!(sim, sim.driver.integrator)
            _apply_torques!()

            run_sim(sim; steps=steps, dt=dt, save_data=dynamic_data_save,
                    save_m_every=dynamic_m_every, call_back=call_back)
        end

        if vtk_saving
            !isdir("vtks") && mkdir("vtks")
            save_vtk(sim, @sprintf("vtks/%s.vts", name))
        end

        return sim
    end

    # Now we need to deal with the case that N > 0
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

        if _is_relax(task_)
            # honor driver_s for relax stages as well (I-05)
            if lowercase(string(driver_)) != lowercase(sim.driver_name)
                set_driver(sim; driver=driver_, integrator=integ, driver_kw...)
            end
            _apply_torques!()

            relax(sim; max_steps=max_steps, stopping_dmdt=stopping_dmdt,
                  save_data_every=relax_data_every, save_m_every=relax_m_every,
                  using_time_factor=using_time_factor)

        elseif _is_dynamics(task_)
            if lowercase(string(driver_)) == "sd"
                # SD cannot run dynamics; fall back to LLG (I-03: re-apply params)
                @warn "driver=\"SD\" cannot run a Dynamics stage; switched to \"LLG\"."
                set_driver(sim; driver="LLG", integrator=integ, driver_kw...)
            elseif lowercase(string(driver_)) != lowercase(sim.driver_name)
                set_driver(sim; driver=driver_, integrator=integ, driver_kw...)
            end
            set_initial_condition!(sim, sim.driver.integrator)
            _apply_torques!()

            run_sim(sim; steps=steps, dt=dt, save_data=dynamic_data_save,
                    save_m_every=dynamic_m_every, call_back=call_back)

        else
            error("Only support two types of task: 'Relax' and 'Dynamics'.")
        end

        if vtk_saving
            !isdir("vtks") && mkdir("vtks")
            vtkname = @sprintf("vtks/%s_%d.vts", name, n)
            save_vtk(sim, vtkname)
        end
    end

    return sim
end

function sim_with(; args...)
    return sim_with(Dict(args))
end
