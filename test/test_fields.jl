using MicroMag
using Test
using DelimitedFiles

function load_field_data()
    filename = joinpath(@__DIR__, "test_fields.txt")
    data = readdlm(filename, ' ', Float64; comments=true)
    #println(size(data))
    N, n = size(data)
    Ms = data[1:Int((N - 1) / 3), 1]
    m0 = data[1:(N - 1), 2]
    demag = data[1:(N - 1), 3]
    exch = data[1:(N - 1), 4]
    dmi = data[1:(N - 1), 5]
    dmi_int = data[1:(N - 1), 6]
    anis = data[1:(N - 1), 7]
    #The last row is enery
    return Ms, m0, demag, exch, dmi, dmi_int, anis, data[N, 3:7]
end

function test_fields(A=1.3e-11, D=4e-3, DI=2e-3, Ku=-3e4, axis=(0, 0, 1))
    mesh = FDMesh(; nx=20, ny=5, nz=3, dx=2.5e-9, dy=2.5e-9, dz=3e-9)
    Ms, m0, fd, fe, fdmi, fdmi2, fan, energy = load_field_data()
    sim = Sim(mesh)
    set_Ms(sim, Ms)
    init_m0(sim, m0; norm=false)
    exch = add_exch(sim, A)
    demag = add_demag(sim)
    dmi = add_dmi(sim, D)
    dmi2 = add_dmi(sim, DI; type="interfacial")
    anis = add_anis(sim, Ku; axis=axis)

    MicroMag.effective_field(sim, sim.spin, 0.0)
    MicroMag.effective_field(sim, sim.spin, 0.0)

    #println("max exch error: ",maximum(abs.(exch.field-fe)))
    #println(maximum(abs.(demag.field-fd)))
    #println(maximum(abs.(dmi.field-fdmi)))
    #println(maximum(abs.(dmi2.field-fdmi2)))
    #println(maximum(abs.(anis.field-fan)))
    @test isapprox(Array(exch.field), fe, atol=1e-5)
    @test isapprox(Array(demag.field), fd, atol=1e-5)
    @test isapprox(Array(dmi.field), fdmi, atol=1e-5)
    @test isapprox(Array(dmi2.field), fdmi2, atol=1e-5)
    @test isapprox(Array(anis.field), fan, atol=1e-5)

    @test (abs((sum(demag.energy) - energy[1]) / energy[1])) < 1e-11
    @test (abs((sum(exch.energy) - energy[2]) / energy[2])) < 1e-15
    @test (abs((sum(dmi.energy) - energy[3]) / energy[3])) < 1e-15
    @test (abs((sum(dmi2.energy) - energy[4]) / energy[4])) < 1e-15
    @test (abs((sum(anis.energy) - energy[5]) / energy[5])) < 1e-15
end

@testset "Test Fields" begin
    set_backend("cpu")
    test_fields()
end

@testset "Test Fields GPU" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        test_fields()
    end
end
