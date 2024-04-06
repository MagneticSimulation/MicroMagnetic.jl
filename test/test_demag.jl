using JuMag
using Test

function test_demag()
    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=2, ny=1, nz=1)
    Ms = 8.6e5
    sim = Sim(mesh)
    set_Ms(sim, Ms)

    init_m0(sim, (1, 0, 0))
    add_demag(sim)

    JuMag.effective_field(sim, sim.spin, 0.0)
    #println(sim.field)
    @test isapprox(Array(sim.field), [-170551.8913984, 0, 0, -170551.8913984, 0, 0])

    mesh = FDMesh(; nx=3, ny=2, nz=1)
    Ms = 8.6e5
    sim = Sim(mesh)
    set_Ms(sim, Ms)
    init_m0(sim, (0.1, 0.2, 1))
    add_demag(sim)

    JuMag.effective_field(sim, sim.spin, 0.0)
    #println(sim.field)
    expected = [-9615.99019074, -39898.78767025, -430282.70478141, -8664.33293854,
                -51323.59117349, -496012.77748287, -27749.99302205, -48965.78908591,
                -430282.70478141, -27749.99302205, -48965.78908591, -430282.70478141,
                -8664.33293854, -51323.59117349, -496012.77748287, -9615.99019074,
                -39898.78767025, -430282.70478141]
   @test isapprox(Array(sim.field), expected)

    mesh = FDMesh(; dx=2, dy=1, dz=1, nx=12, ny=1, nz=1)
    Ms = 8.6e5
    sim = Sim(mesh)
    set_Ms(sim, Ms)

    init_m0(sim, (0.5, 0.6, 0.7))
    add_demag(sim)

    JuMag.effective_field(sim, sim.spin, 0.0)

    expected = [-40715.4448514, -221564.08111452, -258491.42796693, -3885.3038712,
                -243662.16570271, -284272.52665316, -1420.95983327, -245140.77212539,
                -285997.56747963, -785.62732964, -245521.97162757, -286442.30023216,
                -550.58082889, -245662.99952803, -286606.8327827, -464.3653603,
                -245714.72880919, -286667.18361072, -464.3653603, -245714.72880919,
                -286667.18361072, -550.58082889, -245662.99952803, -286606.8327827,
                -785.62732964, -245521.97162757, -286442.30023216, -1420.95983327,
                -245140.77212539, -285997.56747963, -3885.3038712, -243662.16570271,
                -284272.52665316, -40715.4448514, -221564.08111452, -258491.42796693]
    @test isapprox(Array(sim.field), expected)
end

function compute_field_macro_uniform(mesh, Nx, Ny, Nz)
    Ms = 8.6e5
    sim = Sim(mesh)
    sim.Ms[:] .= Ms

    init_m0(sim, (1, 1, 1); norm=false)
    add_demag(sim; Nx=Nx, Ny=Ny, Nz=Nz)

    JuMag.effective_field(sim, sim.spin, 0.0)

    return Array(sim.field) / Ms
end

function compute_field_macro(mesh, Nx, Ny, Nz)
    Ms = 8.6e5
    sim = Sim(mesh)
    sim.Ms[:] .= Ms

    function rand_m(i, j, k, dx, dy, dz)
        I = i % 3
        J = j % 3
        return (0.2, sin(I * 0.9), cos(J * 1.1))
    end

    init_m0(sim, rand_m)
    add_demag(sim; Nx=Nx, Ny=Ny, Nz=Nz)

    JuMag.effective_field(sim, sim.spin, 0.0)

    return Array(sim.field) / Ms
end

function test_field_macro()
    mesh1 = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=1)
    mesh2 = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=3, nz=1)
    mesh3 = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=101, ny=101, nz=1)
    mesh4 = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=33, ny=33, nz=1)

    f1 = compute_field_macro_uniform(mesh1, 500, 0, 0)  #1d case
    #println(f1)
    @test isapprox(f1[1], 0, atol=1e-6)
    @test isapprox(f1[2], -0.5, atol=1e-6)
    @test isapprox(f1[3], -0.5, atol=1e-6)

    f2 = compute_field_macro_uniform(mesh3, 5, 5, 0)  #2d case
    field = reshape(f2, 3, 101, 101, 1)
    f2b = field[:, 51, 51, 1]
    #println(f2b)
    @test isapprox(f2b[1], 0, atol=1e-3)
    @test isapprox(f2b[2], 0, atol=1e-3)
    @test isapprox(f2b[3], -1, atol=1e-3)

    f3 = compute_field_macro_uniform(mesh1, 31, 31, 0) #self-consistency
    f4 = compute_field_macro_uniform(mesh2, 10, 10, 0)
    field = reshape(f4, 3, 3, 3, 1)
    f3b = field[:, 2, 2, 1]
    #println(f3, f3b)
    @test isapprox(f3, f3b, atol=1e-10)

    f5 = compute_field_macro(mesh2, 5, 5, 0) #self-consistency
    f6 = compute_field_macro(mesh4, 0, 0, 0)
    f5 = Array(f5)
    f6 = Array(f6)
    field = reshape(f6, 3, 33, 33, 1)
    f6b = field[1:3, 16:18, 16:18, 1]
    field = reshape(f6b, 27)
    #println(maximum(abs.(field)))
    #println(maximum(abs.(f5 - field)))
    @test maximum(abs.(f5 - field)) < 1e-10
end

@testset "Test Demag CPU" begin
    set_backend("cpu")
    test_demag()
    test_field_macro()
end

@testset "Test Demag CUDA" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        test_demag()
        test_field_macro()
    end
end
