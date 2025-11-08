using MicroMagnetic, Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_sllg(; dt=1e-15)
    mesh = CubicMesh(; nx=3, ny=3, nz=3)

    V = 2.8e-26
    sim = Sim(mesh; driver="LLG", integrator="RungeKutta", name="sllg")
    sim.driver.alpha = 0.1
    sim.driver.gamma = 1.76e11
    sim.driver.integrator.step = dt

    set_mu_s(sim, 1.42e5 * V)

    init_m0_random(sim)

    add_anis(sim, 7.2e5 * V; axis=(0, 0, 1))
    add_thermal_noise(sim, 300.00; scaling=t -> exp(-t / 1e-12))

    relax(sim; max_steps=5000, stopping_dmdt=0, save_data_every=100)

    data, units = read_table("sllg_llg.txt")

    @test abs(data["time"][end] - 5000 * dt) < 1e-15

    eps = MicroMagnetic.Float[] == Float32 ? 1e-7 : 1e-12
    @test abs(data["T_thermal"][end] - 300 * exp(-5000 * dt / 1e-12)) < eps
end

function test_noise(; dt=1e-14)
    mesh = FDMesh(; nx=4, ny=3, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; driver="LLG", integrator="Heun", name="sllg2")
    sim.driver.alpha = 0.1
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.step = dt

    set_Ms(sim, 8e5)

    init_m0_random(sim)

    t1 = add_thermal_noise(sim, 300.0)
    t2 = add_thermal_noise(sim, 100.0; scaling=t -> exp(-t/1e-9), T0=23)

    gaussian_profile = (x, y, z) -> 200 * exp(-(x^2 + y^2)/1e-18)
    t3 = add_thermal_noise(sim, gaussian_profile; T0=50)

    dynamic_temp = (x, y, z, t) -> 300 + 20*sin(2π*t/1e-9 + 0.1*x*1e9)*exp(-t/2e-9)
    t4 = add_thermal_noise(sim, dynamic_temp; name="t4")

    dynamic_temp2 = (x, y, z) -> 300 + 20*sin(0.1*x*1e9)*exp(-0.5)
    t5 = add_thermal_noise(sim, dynamic_temp2; T0=50)

    MicroMagnetic.effective_field(sim, sim.spin, 1e-9)

    @test maximum(t1.temperature) == 300
    @test maximum(t2.temperature) == 100
    @test isapprox(t2.scaling_factor, exp(-1.0))

    t3_temp = Array(t3.temperature)
    @test isapprox(t3_temp[1], 200*exp(-1.5^2-1^2))

    t4_temp = Array(t4.temperature)
    @test isapprox(t4_temp[1], 300 + 20*sin(-0.15)*exp(-0.5))

    @test isapprox(t4_temp, Array(t5.temperature))

    @test t4.name == "t4"
end

@using_gpu()
test_functions("SLLG", test_sllg, test_noise)
