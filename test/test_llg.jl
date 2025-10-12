using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_llg(integrator="DormandPrince")
    tol = integrator == "BS23" ? 1e-7 : 1e-8
    error = integrator == "BS23" ? 2e-6 : 8e-7

    mesh = FDMesh(; nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5

    sim.driver.integrator.tol = tol

    add_zeeman(sim, (0, 0, 1e5))

    init_m0(sim, (1.0, 0, 0))

    for i in 1:10
        run_until(sim, 1e-10 * i)
    end

    #println(sim.driver.integrator.t, " ", integrator)

    ts = Array([1e-9])
    mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)

    m = Array(sim.spin)

    @test abs(mx[1] - m[1]) < error
    @test abs(my[1] - m[2]) < error
    @test abs(mz[1] - m[3]) < error
end

function test_llg_rk(integrator="RungeKutta")
    #Test mesh
    mesh = FDMesh(; nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)
    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.step = 1e-13

    add_zeeman(sim, (0, 0, 1e5))
    init_m0(sim, (1.0, 0, 0))

    for i in 1:500
        MicroMagnetic.run_step(sim, sim.driver)
    end

    #println(sim.spin[1]," ",sim.spin[2]," ",sim.spin[3])
    ts = Array([5e-11])
    mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)

    m = Array(sim.spin)

    @test abs(mx[1] - m[1]) < 8e-6
    @test abs(my[1] - m[2]) < 8e-6
    @test abs(mz[1] - m[3]) < 8e-6
end

function test_llg_dp()
    test_llg("DormandPrince")
    return test_llg("DormandPrinceCayley")
end

function test_llg_bs()
    return test_llg("BS23")
end

function test_llg_ck()
    return test_llg("CashKarp54")
end

function test_llg_fh()
    return test_llg("Fehlberg54")
end

function test_llg_fixed_dt()
    for integrator in ["Heun", "RungeKutta", "RungeKuttaCayley"]
        test_llg_rk(integrator)
    end
end

@using_gpu()
test_functions("LLG DormandPrince", test_llg_dp)
test_functions("LLG BS23", test_llg_bs)
test_functions("LLG CashKarp54", test_llg_ck)
test_functions("LLG Fehlberg54", test_llg_fh)
test_functions("LLG Fixed dt", test_llg_fixed_dt)
