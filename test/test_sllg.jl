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

#@using_gpu()
test_functions("SLLG", test_sllg)
