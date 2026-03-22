using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function init_m_fun(x, y, z)
    r = sqrt(x * x + y * y + z * z)
    return (0, sin(r), cos(r))
end

function test_interpolate()
    filepath = joinpath(@__DIR__, "meshes/cylinder.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5)

    init_m0(sim, (0,0,1))

    spin = Array(sim.spin)
    point = [0.1, 0.2, 2.0]
    points = hcat([0.5, 0.4, 1.1], [1.6, 0.5, 2.4])

    f1 = interpolate_field(mesh, spin, point)
    f2 = interpolate_field(mesh, spin, points)
    println(f1) 
    println(f2)
    @test f1 ≈ [0,0,1]
    @test f2 ≈ reshape([0,0,1, 0,0,1], (3, 2))
end

test_functions("Test interpolate (FE)", test_interpolate, precisions=[Float64])