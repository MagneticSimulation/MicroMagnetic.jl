using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function time_fun(t)
    return (cos(1e6 * t), cos(1e9 * t), 1)
end

function test_zeeman()
    mesh = FDMesh(; dx=2e-9, nx=3, ny=2, nz=1, pbc="x")

    sim = Sim(mesh; name="test_zeeman")
    sim.driver.alpha = 0.01
    sim.driver.gamma = 2.21e5

    set_Ms(sim, 8.6e5)
    init_m0(sim, (1, 1, 1); norm=false)

    z1 = add_zeeman(sim, (1, 2, 2e3))
    z2 = add_zeeman(sim, (1, 2, 2e3), time_fun)

    T = MicroMagnetic.Float[]
    
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    f1 = Array(z1.field)
    f2 = Array(z2.field)
    @test f1[1] == T(1.0)
    @test f1[2] == T(2.0)
    @test f1[3] == T(2e3)
    @test isapprox(f1, f2)

    @test abs(sum(Array(z1.energy))) > 0
    @test sum(Array(z1.energy)) == sum(Array(z2.energy))


    update_zeeman(sim, (3, 4.9, -10.01))
    f1 = Array(z1.field)
    @test f1[1] == T(3.0)
    @test f1[2] == T(4.9)
    @test f1[3] == T(-10.01)

    z2 = add_zeeman(sim, (1e3, 1e4, 1e5), time_fun)

    MicroMagnetic.effective_field(sim, sim.spin, 1.23e-11)
    f2 = Array(z2.field)
    @test f2[1] == T(1e3 * cos(1.23e-11*1e6))
    @test f2[2] == T(1e4 * cos(1.23e-2))
    @test f2[3] == T(1e5)
    return nothing
end

@using_gpu()
test_functions("Zeeman", test_zeeman)

