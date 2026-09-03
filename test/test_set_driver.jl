using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

# Small sim with a zeeman field so every driver actually produces a nonzero torque;
# data saving is disabled per-call, so the tests leave no artifacts.
function make_sim(; name, driver="LLG")
    mesh = FDMesh(; nx=4, ny=1, nz=1, dx=2e-9)
    sim = Sim(mesh; name=name, driver=driver)
    set_Ms(sim, 8.6e5)
    add_zeeman(sim, (0, 0, 1e5))
    init_m0(sim, (1, 0, 0))
    return sim
end

function relax_quiet(sim; max_steps=2)
    relax(sim; max_steps=max_steps, save_data_every=-1, save_m_every=-1)
end

function check_unit(sim; atol)
    m = reshape(Array(sim.spin), 3, sim.n_total)
    @test all(isfinite, m)
    @test maximum(abs.(vec(sum(abs2, m; dims=1)) .- 1)) < atol
end

# A1 regression: switching *to* SD must not touch the (nonexistent) integrator field
function test_set_driver_to_sd()
    sim = make_sim(name="test_set_driver_sd")
    set_driver(sim; driver="SD")
    @test sim.driver_name == "SD"
    @test sim.driver isa MicroMagnetic.SD
    relax_quiet(sim)
    check_unit(sim; atol=1e-3)
    @test Array(sim.spin)[3] > 0  # tilted toward +z by the zeeman field
end

# Every switch direction, including both directions that enter a driver without
# an integrator field (SD) and the InertialLLG segment via run_until (relax has
# no run_step for InertialLLG, see SET_DRIVER_FIX.md §5.3)
function test_set_driver_roundtrip()
    sim = make_sim(name="test_set_driver_rt", driver="SD")
    relax_quiet(sim)
    set_driver(sim; driver="LLG")
    relax_quiet(sim)
    check_unit(sim; atol=1e-5)
    set_driver(sim; driver="SD")
    relax_quiet(sim)
    check_unit(sim; atol=1e-3)
    set_driver(sim; driver="InertialLLG")
    run_until(sim, 1e-11; save_data=false)
    @test sim.driver_name == "InertialLLG"
    check_unit(sim; atol=1e-3)
end

# EmptyDriver has no integrator field either -- same crash family as SD
function test_set_driver_to_none()
    sim = make_sim(name="test_set_driver_none")
    set_driver(sim; driver="None")
    @test sim.driver_name == "None"
    @test sim.driver isa MicroMagnetic.EmptyDriver
end

# kwargs survive a driver switch and land on the freshly created driver
function test_set_driver_kwargs_on_switch()
    sim = make_sim(name="test_set_driver_kw", driver="SD")
    set_driver(sim; driver="LLG", alpha=0.05, gamma=2.15e5, tol=1e-9)
    @test sim.driver_name == "LLG"
    @test isapprox(first(sim.driver.alpha), 0.05; rtol=1e-5)
    @test isapprox(sim.driver.gamma, 2.15e5; rtol=1e-5)
    @test isapprox(sim.driver.integrator.tol, 1e-9; rtol=1e-5)
end

# Same-name set_driver is a no-op for driver and saver, but kwargs still apply
function test_set_driver_repeat()
    sim = make_sim(name="test_set_driver_repeat")
    d0 = sim.driver
    saver0 = sim.saver
    @test set_driver(sim; driver="LLG") === nothing
    @test sim.driver === d0
    @test sim.saver === saver0
    set_driver(sim; driver="LLG", alpha=0.02, tol=1e-9)
    @test sim.driver === d0
    @test isapprox(first(sim.driver.alpha), 0.02; rtol=1e-5)
    @test isapprox(sim.driver.integrator.tol, 1e-9; rtol=1e-5)
end

# Saver file naming (sd/llg/illg tags), reset on switch, time column per driver
function test_set_driver_saver()
    sim = make_sim(name="test_set_driver_saver", driver="SD")
    @test sim.saver.name == "test_set_driver_saver_sd.txt"
    @test !any(item -> item.name == "time", sim.saver.items)

    set_driver(sim; driver="LLG")
    @test sim.saver.name == "test_set_driver_saver_llg.txt"
    @test any(item -> item.name == "time", sim.saver.items)
    @test sim.saver.t == 0
    @test sim.saver.nsteps == 0
    @test sim.saver.header_saved == false

    set_driver(sim; driver="InertialLLG")
    @test sim.saver.name == "test_set_driver_saver_illg.txt"

    mesh = FDMesh(; nx=4, ny=1, nz=1, dx=2e-9)
    sim2 = Sim(mesh; name="test_set_driver_saver2", driver="InertialLLG")
    @test sim2.saver.name == "test_set_driver_saver2_illg.txt"
    @test any(item -> item.name == "time", sim2.saver.items)
end

# Direct set_driver_arguments on LLG: scalar/array/function alpha, consumption
# contract (applied keys are removed, unsupported keys stay)
function test_set_driver_args_llg()
    sim = make_sim(name="test_set_driver_args_llg")
    args = Dict(:alpha => 0.03, :gamma => 2.15e5, :tol => 1e-8, :unknown => 1)
    MicroMagnetic.set_driver_arguments(sim, args)
    @test !haskey(args, :alpha)
    @test !haskey(args, :gamma)
    @test !haskey(args, :tol)
    @test haskey(args, :unknown)
    @test isapprox(first(sim.driver.alpha), 0.03; rtol=1e-5)
    @test isapprox(sim.driver.gamma, 2.15e5; rtol=1e-5)
    @test isapprox(sim.driver.integrator.tol, 1e-8; rtol=1e-5)

    MicroMagnetic.set_driver_arguments(sim, Dict(:alpha => fill(0.04, sim.n_total)))
    alpha = Array(sim.driver.alpha)  # dense alpha lives on the backend under CUDA
    @test isapprox(alpha[1], 0.04; rtol=1e-5)
    @test eltype(sim.driver.alpha) == eltype(sim.driver.gamma)

    f(i, j, k, dx, dy, dz) = i == 1 ? 0.01 : 0.02
    MicroMagnetic.set_driver_arguments(sim, Dict(:alpha => f))
    alpha = Array(sim.driver.alpha)
    @test isapprox(alpha[1], 0.01; rtol=1e-5)
    @test isapprox(alpha[4], 0.02; rtol=1e-5)
end

# SD consumes nothing (driver type dispatch; no string matching involved)
function test_set_driver_args_sd()
    sim = make_sim(name="test_set_driver_args_sd", driver="SD")
    args = Dict(:alpha => 0.03, :gamma => 2.15e5, :tol => 1e-8)
    MicroMagnetic.set_driver_arguments(sim, args)
    @test haskey(args, :alpha)
    @test haskey(args, :gamma)
    @test haskey(args, :tol)
end

# A2 regression: tol with an integrator that has no tol field must not throw
function test_set_driver_gpsm_tol()
    sim = make_sim(name="test_set_driver_gpsm", driver="SD")
    set_driver(sim; driver="LLG", integrator="GPSM", tol=1e-8)
    @test sim.driver_name == "LLG"
    @test sim.driver.integrator isa MicroMagnetic.GPSM
    args = Dict(:tol => 1e-8)
    MicroMagnetic.set_driver_arguments(sim, args)
    @test haskey(args, :tol)  # GPSM has no tol field: not consumed
end

# The deprecated "SpatialLLG" alias normalizes to LLG, so it is a no-op here
function test_set_driver_spatialllg()
    sim = make_sim(name="test_set_driver_spatialllg")
    d0 = sim.driver
    set_driver(sim; driver="SpatialLLG")
    @test sim.driver === d0
    @test sim.driver_name == "LLG"
end

# R1 regression: InertialLLG accepts alpha/gamma/tol now instead of dropping them
function test_set_driver_illg_kwargs()
    sim = make_sim(name="test_set_driver_illg_kw")
    set_driver(sim; driver="InertialLLG", alpha=0.2, gamma=2.1e5, tol=1e-8)
    @test sim.driver isa MicroMagnetic.InertialLLG
    @test isapprox(sim.driver.alpha, 0.2; rtol=1e-5)
    @test isapprox(sim.driver.gamma, 2.1e5; rtol=1e-5)
    @test isapprox(sim.driver.integrator.tol, 1e-8; rtol=1e-5)
end

# Public API path: driver_s stages that switch *to* SD (P1 before the fix)
function test_sim_with_switch_to_sd()
    mesh = FDMesh(; nx=2, ny=2, nz=1)
    sim = sim_with((name="test_switch_sd", task_s=["relax", "relax"],
                   driver_s=["LLG", "SD"], mesh=mesh, m0=(1, 0, 0), Ms=8e5, A=1e-11,
                   H=(0, 0, 1e5), alpha=0.05, max_steps=2, stopping_dmdt=1e-5,
                   relax_data_every=-1, relax_m_every=-1))
    @test sim.driver_name == "SD"
    @test sim.driver isa MicroMagnetic.SD
end

# SD cannot run a Dynamics stage: coerced to LLG with the driver parameters
# re-applied (I-03), and the coercion is announced with a @warn (R7)
function test_sim_with_dyn_sd_coercion()
    mesh = FDMesh(; nx=2, ny=2, nz=1)
    args = (name="test_dyn_sd", task="Dynamics", driver="SD", mesh=mesh, m0=(1, 0, 0),
            Ms=8e5, H=(0, 0, 1e5), alpha=0.05, steps=1, dt=1e-12,
            dynamic_data_save=false, dynamic_m_every=-1)
    # match_mode=:any: sim_with also emits @info records, and the default :all
    # mode would require the captured logs to be exactly the pattern
    sim = @test_logs match_mode=:any (:warn, "driver=\"SD\" cannot run a Dynamics stage; switched to \"LLG\".") sim_with(args)
    @test sim.driver_name == "LLG"
    @test sim.driver isa MicroMagnetic.LLG
    @test isapprox(first(sim.driver.alpha), 0.05; rtol=1e-5)
end

@using_gpu()
test_functions("set_driver_to_sd", test_set_driver_to_sd)
test_functions("set_driver_roundtrip", test_set_driver_roundtrip)
test_functions("set_driver_to_none", test_set_driver_to_none)
test_functions("set_driver_kwargs_on_switch", test_set_driver_kwargs_on_switch)
test_functions("set_driver_repeat", test_set_driver_repeat)
test_functions("set_driver_saver", test_set_driver_saver)
test_functions("set_driver_args_llg", test_set_driver_args_llg)
test_functions("set_driver_args_sd", test_set_driver_args_sd)
test_functions("set_driver_gpsm_tol", test_set_driver_gpsm_tol)
test_functions("set_driver_spatialllg", test_set_driver_spatialllg)
test_functions("set_driver_illg_kwargs", test_set_driver_illg_kwargs)
test_functions("sim_with_switch_to_sd", test_sim_with_switch_to_sd)
test_functions("sim_with_dyn_sd_coercion", test_sim_with_dyn_sd_coercion)
