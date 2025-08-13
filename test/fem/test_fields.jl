using MicroMagnetic
using Test
using DelimitedFiles

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function test_zeeman()
    filepath = joinpath(@__DIR__, "meshes/octa.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8.6e5)

    z = add_zeeman(sim, (1, 2, 1e5))

    MicroMagnetic.effective_field(sim, sim.spin, 1.23e-11)

    @test z.field[1] == 1.0
    @test z.field[2] == 2.0
    @test z.field[3] == 1e5
end

function init_m_fun(x, y, z)
    r = sqrt(x * x + y * y + z * z)
    return (0, sin(r), cos(r))
end

function test_init_m_function()
    filepath = joinpath(@__DIR__, "meshes/cylinder.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)
    N = length(sim.spin)
    println("The length of spin: ", N)

    filepath = joinpath(@__DIR__, "fields/m.txt")
    m0 = readdlm(filepath; header=true)[1]
    mxyz = reshape(transpose(m0[:, 4:6]), N, 1)
    #println("max diff: ", maximum(abs.(sim.spin - mxyz)))
    eps = 1e-6
    @test maximum(abs.(sim.spin - mxyz)) < eps
    #print(mxyz[1:50])
end

function test_zeeman_field()
    filepath = joinpath(@__DIR__, "meshes/cylinder.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_zeeman(sim, (0, 0, 1e5))

    MicroMagnetic.effective_field(z, sim, sim.spin, 0.0)

    expected_energy = 1.81151949357e-21
    println("zeeman energy: ", sum(z.energy))
    @test abs(sum(z.energy) - expected_energy) / expected_energy < 1e-10
end

function test_anis_field()
    filepath = joinpath(@__DIR__, "meshes/cylinder.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_anis(sim, 5.2e5; axis=(0, 0.6, 0.8))

    MicroMagnetic.effective_field(z, sim, sim.spin, 0.0)

    N = 3 * sim.n_total

    filepath = joinpath(@__DIR__, "fields/anis.txt")
    f0 = readdlm(filepath; header=true)[1]
    field = reshape(transpose(f0[:, 4:6]), N, 1)
    eps = 1e-6
    @test maximum(abs.(Array(z.field) - field)) < eps

    #expected_energy = 9.60333028541e-21

    #println("anis energy: ",sum(z.energy))
    #@test (sum(z.energy)-expected_energy)/expected_energy<1e-10

end

function test_exchange_field()
    filepath = joinpath(@__DIR__, "meshes/cylinder.mesh")
    mesh = FEMesh(filepath)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_exch(sim, 1.3e-11)

    MicroMagnetic.effective_field(z, sim, sim.spin, 0.0)

    N = 3 * sim.n_total

    filepath = joinpath(@__DIR__, "fields/exch.txt")
    f0 = readdlm(filepath; header=true)[1]
    field = reshape(transpose(f0[:, 4:6]), N, 1)
    eps = 5e-5
    @test maximum(abs.(Array(z.field) - field)) < eps

    expected_energy = 3.18671455452e-19

    println("exch energy: ", sum(z.energy))
    @test abs(sum(z.energy) - expected_energy) / expected_energy < 1e-10
end



test_zeeman()
test_init_m_function()

@using_gpu()
test_functions("Test zeeman (FE)", test_zeeman_field, precisions=[Float64])
test_functions("Test anisotropy (FE)", test_anis_field, precisions=[Float64])
test_functions("Test exchange (FE)", test_exchange_field, precisions=[Float64])