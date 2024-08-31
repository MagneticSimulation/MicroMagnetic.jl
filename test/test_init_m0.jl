using MicroMagnetic
using Test
using LinearAlgebra

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function init_fun3(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end

function init_fun6(i, j, k, dx, dy, dz)
    x = i - 50.5
    y = j - 50.5
    r = (x^2 + y^2)^0.5
    if r < 20
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end


function test_init_m0()
    mesh = FDMesh(; nx=100, ny=100, nz=2)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)

    init_m0(sim, init_fun3)

    m1 = Array(sim.spin)

    init_m0(sim, init_fun6)
    m2 = Array(sim.spin)

    @test isapprox(m1, m2)

end


@using_gpu()
test_functions("init_m0", test_init_m0)