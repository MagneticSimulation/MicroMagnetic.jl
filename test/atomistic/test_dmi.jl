using MicroMagnetic
using Test
using DelimitedFiles

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function spatial_bulk_DMI(i, j, k)
    Dx = 0.3meV
    Dy = 0.3meV
    Dz = 0.3meV
    return [(Dx, 0, 0), (0, Dy, 0), (0, 0, Dz)]
end

function spatial_m0(i, j, k, dx, dy, dz)
    return (sin(i / 10), cos(i / 10+ 5*j/10), 0.3)
end

function test_atomistic_dmi()
    mesh = CubicMesh(; nx=12, ny=4, nz=1, dx=0.5e-9, dy=0.5e-9, dz=0.5e-9, pbc="y")

    sim = Sim(mesh)

    set_mu_s(sim, mu_B)

    init_m0(sim, spatial_m0)

    dmi1 = add_dmi(sim, 0.3meV)
    dmi2 = add_dmi(sim, spatial_bulk_DMI)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    @test isapprox(Array(dmi1.field), Array(dmi2.field), atol=1e-7)
end

#test_atomistic_dmi()

@using_gpu()
test_functions("Atomistic DMI", test_atomistic_dmi, precisions=[Float64])
