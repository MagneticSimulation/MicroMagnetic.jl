using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function init_saving(mesh, m)
    sim = Sim(mesh, name = "saving")
    set_Ms(sim, 8.0e5)
    init_m0(sim, m)
    save_ovf(sim, "Float64", type = Float64)
    save_ovf(sim, "Float32", type = Float32)
    save_ovf(sim, "String", type = String)
end

function test_reading(fname, mesh, m)
    ovf = read_ovf(fname)
    @test length(ovf.data) == 6
    @test maximum(abs.(ovf.data-m)) < 1e-6
end

function test_sim_reading(fname, mesh, m)
    sim = Sim(mesh, name = "reading")
    set_Ms(sim, 8.0e5)
    read_ovf(sim, fname)
    @test isapprox(Array(sim.spin), m, atol=1e-6)
end

function test_copying(fname, mesh, m)
    ovf = read_ovf(fname)
    save_ovf(ovf, fname*"_copy_float64")

    ovf = read_ovf(fname*"_copy_float64")
    @test length(ovf.data) == 6
    @test maximum(abs.(ovf.data-m)) < 1e-6
    rm(fname*"_copy_float64.ovf")

    ovf = read_ovf(fname)
    save_ovf(ovf, fname*"_copy_string", type = String)
    ovf = read_ovf(fname*"_copy_string")
    @test length(ovf.data) == 6
    @test maximum(abs.(ovf.data-m)) < 1e-6
    rm(fname*"_copy_string.ovf")
end

function test_mag2ovf(m::Array, nx::Int, ny::Int, nz::Int)
    mag2ovf(m, nx, ny, nz, fname = "test")
    ovf = read_ovf("test")
    @test(isapprox(ovf.data, m, atol=1e-6))
    @test(ovf.xnodes == nx)
    @test(ovf.ynodes == ny)
    @test(ovf.znodes == nz)
    rm("test.ovf")
end

function test_ovf()
    m = [0.6, 0.8, 0, 0.6, 0.8, 0]
    mesh = FDMesh(nx = 2, ny = 1, nz = 1, dx = 1e-9, dy = 1e-9, dz = 1e-9)
    init_saving(mesh, m)
    for fname in ["Float64", "Float32", "String"]
        test_reading(fname, mesh, m)
        test_sim_reading(fname, mesh, m)
        test_copying(fname, mesh, m)
    end

    test_mag2ovf(m, 2, 1, 1)
end

@using_gpu()
test_functions("OVF", test_ovf)