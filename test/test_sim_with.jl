using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function analytical_llg(alpha::Float64, gamma::Float64, H0::Float64, t::Float64)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * t

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return [mx, my, mz]
end

function test_sim_with()
    mesh = FDMesh(; nx=2, ny=2, nz=2)
    args = (task="Relax", mesh=mesh, m0=(1, 2, 3), Ms=8e5, A=1e-11, H=(0, 0, 1e5), 
            relax_m_every=-1,
            name="test_relax", stopping_dmdt=0.01)

    sim = sim_with(args)
    m = Array(sim.spin)
    @test isapprox(m[1:3], [0, 0, 1])

    mesh = FDMesh(; nx=1, ny=1, nz=1)
    analytical = analytical_llg(0.05, 2.21e5, 1e5, 3e-10)
    args = (task="Dynamics", mesh=mesh, m0=(1, 0, 0), Ms=8e5, H=(0, 0, 1e5), gamma=2.21e5,
            alpha=0.05, steps=10, dt=3e-11, tol=1e-8)

    sim = sim_with(args)
    m = Array(sim.spin)
    @test isapprox(m[1:3], analytical)

    args = (task_s=["relax", "dynamics"], mesh=mesh, m0=(1, 2, 3), Ms=8e5, A=1e-11, steps=10,
            alpha=0.05, tol=1e-8, dt=3e-11, stopping_dmdt=0.001,
            H_s=[(1e5, 0, 0), (0, 0, 1e5)])

    sim = sim_with(args)
    m = Array(sim.spin)
    println(m[1:3], analytical)
    @test isapprox(m[1:3], analytical)

    args = (name="test_read", task="relax", mesh=mesh, m0=(1, 2, 3), Ms=8e5, A=1e-11, stopping_dmdt=0.01,
            relax_m_every =-1,
            H_s=[(1e5, 0, 0), (0, 0, -1e5)])

    sim = sim_with(args)
    m = Array(sim.spin)
    println(m[1:3], [0, 0, -1])
    @test isapprox(m[1:3], [0, 0, -1])
    return nothing
end

function test_stt_kwarg()
    mesh = FDMesh(; nx=2, ny=2, nz=1)
    # I-10: stt kwarg is forwarded to add_stt and applied for LLG-family stages
    args = (name="test_stt_kwarg", task_s=["relax", "dynamics"], driver_s=["SD", "LLG"],
            mesh=mesh, m0=(1, 0, 0), Ms=8e5, H=(0, 0, 1e5), alpha=0.05,
            stt=(model=:zhang_li, b=-10.0, J=(1, 0, 0), xi=0.05),
            max_steps=2, stopping_dmdt=1e-5, steps=1, dt=1e-12,
            relax_m_every=-1, dynamic_m_every=-1)
    sim = sim_with(args)
    @test sim.driver_name == "LLG"
    @test count(i -> i isa MicroMagnetic.ZhangLiTorque, sim.interactions) == 1

    # pure-SD run: stt is never applied to the energy minimization
    sim = sim_with((name="test_stt_sd", task="Relax", mesh=mesh, m0=(1, 0, 0), Ms=8e5,
                    H=(0, 0, 1e5), stt=(model=:zhang_li, b=-10.0, J=(1, 0, 0), xi=0.05),
                    max_steps=1, stopping_dmdt=0.5, relax_m_every=-1))
    @test isempty(filter(i -> i isa MicroMagnetic.ZhangLiTorque, sim.interactions))

    # driver="LLG" at creation: stt applied immediately
    sim = sim_with((name="test_stt_llg", task="Relax", driver="LLG", mesh=mesh,
                    m0=(1, 0, 0), Ms=8e5, H=(0, 0, 1e5), alpha=0.05,
                    stt=(model=:zhang_li, b=-10.0, J=(1, 0, 0), xi=0.05),
                    max_steps=1, stopping_dmdt=0.5, relax_m_every=-1))
    @test count(i -> i isa MicroMagnetic.ZhangLiTorque, sim.interactions) == 1
    return nothing
end

function test_read_table()
        data, units = read_table("test_read_sd.txt")
        @test data["m_x"][end] < 1e-15
        @test data["m_y"][end] < 1e-15
        @test data["m_z"][end]+1 < 1e-15
end

function test_sim_with_validation()
    mesh = FDMesh(; nx=1, ny=1, nz=1)
    # I-04: an invalid task value errors instead of silently doing nothing
    @test_throws ErrorException sim_with((task="minimize", mesh=mesh, m0=(1, 2, 3),
                                          Ms=8e5, name="test_task"))
    # I-09: unknown keys error before the simulation starts
    err = try
        sim_with((mesh=mesh, alpah=0.05, m0=(1, 2, 3), Ms=8e5, name="test_keys",
                  task="Relax", stopping_dmdt=0.5, max_steps=1))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("alpha", err.msg)
    # I-07: sweeping over keys without per-stage support errors
    @test_throws ErrorException sim_with((mesh=mesh, alpha_s=[0.01, 0.1], m0=(1, 2, 3),
                                          Ms=8e5, task="Relax", max_steps=1,
                                          stopping_dmdt=0.5, name="test_sweep_key"))
    # I-08: the caller's Dict is not modified
    args = Dict(:mesh => mesh, :Ms => 8e5, :m0 => (1, 2, 3), :name => "test_dict",
                :A => 1e-11, :stopping_dmdt => 0.5, :max_steps => 2, :relax_m_every => -1)
    sim_with(args)
    @test args[:Ms] == 8e5
    @test haskey(args, :A)
    @test args[:stopping_dmdt] == 0.5
    return nothing
end

function test_sim_with_sweep_driver()
    mesh = FDMesh(; nx=1, ny=1, nz=1)
    # I-03: driver parameters survive driver switches between stages
    args = (name="test_sweep_driver", task_s=["dynamics", "dynamics"],
            driver_s=["LLG", "LLG_STT"], mesh=mesh, m0=(1, 0, 0), Ms=8e5, H=(0, 0, 1e5),
            alpha=0.05, gamma=2.21e5, steps=1, dt=1e-12, dynamic_m_every=-1)
    sim = sim_with(args)
    @test sim.driver_name == "LLG_STT"
    @test isapprox(sim.driver.alpha, 0.05; rtol=1e-6)
    @test isapprox(sim.driver.gamma, 2.21e5; rtol=1e-6)

    # I-05: relax stages honor driver_s as well
    args = (name="test_sweep_relax", task_s=["relax", "relax"], driver_s=["SD", "LLG"],
            mesh=mesh, m0=(1, 0, 0), Ms=8e5, H=(0, 0, 1e5), alpha=0.05,
            max_steps=2, stopping_dmdt=1e-5, relax_m_every=-1)
    sim = sim_with(args)
    @test sim.driver_name == "LLG"
    @test isapprox(sim.driver.alpha, 0.05; rtol=1e-6)
    return nothing
end

function test_save_at_end()
    mesh = FDMesh(; nx=2, ny=2, nz=1)
    # I-06: relax reaching max_steps without convergence still saves the final state
    args = (name="test_relax_end", task="Relax", mesh=mesh, m0=(1, 0, 0), Ms=8e5,
            A=1e-11, H=(0, 0, 1e5), alpha=0.05, max_steps=3, stopping_dmdt=1e-12,
            relax_data_every=0, relax_m_every=0)
    sim_with(args)
    @test isfile(joinpath("test_relax_end_SD", "m_00000003.ovf"))
    data, units = read_table("test_relax_end_sd.txt")
    @test length(data["m_z"]) == 1

    # I-06: run_sim with save_m_every=0 saves the final state
    sim_with((name="test_dyn_end", task="Dynamics", mesh=mesh, m0=(1, 0, 0), Ms=8e5,
              H=(0, 0, 1e5), alpha=0.05, steps=2, dt=1e-12, dynamic_m_every=0))
    @test isfile(joinpath("test_dyn_end_LLG", "m_00000002.ovf"))
    return nothing
end


@using_gpu()
test_functions("sim_with", test_sim_with)
test_functions("sim_with_validation", test_sim_with_validation)
test_functions("sim_with_sweep_driver", test_sim_with_sweep_driver)
test_functions("save_at_end", test_save_at_end)
test_functions("stt_kwarg", test_stt_kwarg)
test_functions("read_table", test_read_table)