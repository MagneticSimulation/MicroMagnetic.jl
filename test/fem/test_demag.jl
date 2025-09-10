using MicroMagnetic
using Test
using DelimitedFiles
if Base.find_package("CUDSS") !== nothing
    using CUDSS
else
    @warn("CUDSS is needed for FE demag when using CUDA.")
end

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function potential(x, y, z)
    R = sqrt(x * x + y * y + z * z)
    R0 = 10
    M = 8e5
    if R > R0
        return z*M*(R0/R)^3/3.0
    end
    return z*M/3.0
end

function test_demag_field_sphere(; method="bem")
    filepath = joinpath(@__DIR__, "meshes/sphere_air.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5, region_id=1)

    init_m0(sim, (0,0,1))

    demag = add_demag(sim)

    v = zeros(sim.n_total)
    MicroMagnetic.init_scalar_nodes!(v, mesh, potential)
    
    MicroMagnetic.effective_field(sim, sim.spin)
    N = 3 * sim.n_total

    phi = Array(demag.phi1)

    id = argmax(abs.(v .- phi))

    error = abs((v[id] - phi[id])/v[id])

    println(error)
    #save_vtk(sim, "m_demag.vtu", fields=["demag"])

    @test error < 0.04
end

function test_exch_sphere_air()
    filepath = joinpath(@__DIR__, "meshes/sphere_air.mesh")
    mesh = FEMesh(filepath, unit_length=1e-9)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5, region_id=1)

    init_m0(sim, (0,0,1))

    f = add_exch(sim, 1e-13)

    MicroMagnetic.effective_field(sim, sim.spin)

    @test !any(isnan, f.field)
    @test maximum(abs.(f.field)) < 1e-8
end

@using_gpu()
test_functions("Test Demag (FE)", test_demag_field_sphere, precisions=[Float64])
test_functions("Test Exch Air (FE)", test_exch_sphere_air, precisions=[Float64])