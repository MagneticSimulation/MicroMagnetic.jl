using MicroMagnetic
using Test
using DelimitedFiles

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function spatial_Js(i, j, k)
    J = 50 * k_B
    return [J, J, J]
end
function spatial_m0(i, j, k, dx, dy, dz)
    return (sin(i / 10), cos(i / 10), 0.3)
end

function test_atomistic_exch()
    mesh = CubicMesh(; nx=12, ny=4, nz=1, dx=0.5e-9, dy=0.5e-9, dz=0.5e-9, pbc="y")

    sim = Sim(mesh)

    set_mu_s(sim, mu_B)

    init_m0(sim, spatial_m0)

    J = 50 * k_B
    ex1 = add_exch(sim, J)
    ex2 = add_exch(sim, spatial_Js)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    @test isapprox(Array(ex1.field), Array(ex2.field), atol=1e-7)
end

@using_gpu()
test_functions("Atomistic Exch", test_atomistic_exch; precisions=[Float64])