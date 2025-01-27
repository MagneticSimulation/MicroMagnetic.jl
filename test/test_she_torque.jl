using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_she_torque(integrator="DormandPrince")
    #Test mesh
    mesh = FDMesh(; nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.tol = 1e-8

    add_zeeman(sim, (0, 0, 1e5))

    t = add_she_torque(sim, (0.1, 0.2, 0.3), (0.3, 0.4, 0.5), a1=(0.1, 0.2, 0.8), a2=(0.5, -0.9, 0.2), beta=0.135)

    init_m0(sim, (-0.3, 0.4, 0.8), norm=false)
    
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    expected = [-5.619678515837106, 8.710639746606335, -7.754204488687782]
    @test isapprox(Array(t.field)*1e7, expected, rtol=1e-6)
end

@using_gpu()
test_functions("SHETorque", test_she_torque)