using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_sim()
    #Test mesh
    mesh = FDMesh(; dx=1.1e-9, nx=10, ny=2)

    sim = Sim(mesh; driver="LLG")
    set_Ms(sim, 1.0)

    @test sim.n_total == 20

    init_m0(sim, (1.0, 1.0, 0))
    #println(sim.spin)
    spin = Array(sim.spin)

    @test isapprox(spin[1], 0.5 * 2^0.5)
    @test isapprox(spin[2], 0.5 * 2^0.5)

    set_Ms(sim, 8.6e5)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.21e5
    sim.driver.precession = false

    add_exch(sim, 1.3e-11)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    @info("test_sim() passed!")
    #println(sim.spin)
end

@using_gpu()
test_functions("Sim", test_sim)
