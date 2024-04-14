using JuMag
using Test
using DelimitedFiles

function load_field_data()
    filename = joinpath(@__DIR__, "test_fields_atomistic.txt")
    data = readdlm(filename, ' ', Float64; comments=true)
    N, n = size(data)
    Ms = data[1:Int((N - 1) / 3), 1]
    m0 = data[1:(N - 1), 2]
    demag = data[1:(N - 1), 3]
    exch = data[1:(N - 1), 4]
    dmi = data[1:(N - 1), 5]
    anis = data[1:(N - 1), 6]
    #The last row is energy
    return Ms, m0, demag, exch, dmi, anis, data[N, 3:6]
end

function test_atomistic_fields(A=1.3e-11, D=4e-3, DI=2e-3, Ku=-3e4, axis=(0, 0, 1))
    mesh = CubicMesh(; nx=39, ny=11, nz=3, dx=0.5e-9, dy=0.5e-9, dz=0.5e-9)

    mu_s, m0, fd, fe, fdmi, fan, energy = load_field_data()
    sim = Sim(mesh)

    set_mu_s(sim, mu_s)

    init_m0(sim, m0; norm=false)

    J = 50 * k_B
    exch = add_exch(sim, J)
    demag = add_demag(sim)
    dmi = add_dmi(sim, 0.09 * J)
    anis = add_anis(sim, 5e-3 * J; axis=(0, 0, 1))

    JuMag.effective_field(sim, sim.spin, 0.0)

    @test isapprox(Array(exch.field), fe, atol=1e-8)
    @test isapprox(Array(demag.field), fd, atol=1e-8)
    @test isapprox(Array(dmi.field), fdmi, atol=1e-8)
    @test isapprox(Array(anis.field), fan, atol=1e-8)

    delta_anis = 5e-3 * J * sum(mu_s .> 0)
    @test (abs((sum(demag.energy) - energy[1]) / energy[1])) < 1e-11
    @test (abs((sum(exch.energy) - energy[2]) / energy[2])) < 1e-14
    @test (abs((sum(dmi.energy) - energy[3]) / energy[3])) < 1e-14
    @test (abs((sum(anis.energy) - delta_anis - energy[4]) / energy[4])) < 1e-14
end

@testset "Test Atomistic Fields CPU" begin
    set_backend("cpu")
    test_atomistic_fields()
end

@testset "Test Atomistic Fields CUDA" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        test_atomistic_fields()
    end
end
