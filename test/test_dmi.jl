using NuMag
using Test

function test_bulk_dmi()
    function m0_fun(i, j, k, dx, dy, dz)
        if i == 1
            return (1, 2, 3)
        elseif i == 2
            return (4, 5, 6)
        end
        return (1, 0, 0)
    end

    mesh = FDMesh(; dx=2e-9, nx=2, ny=1, nz=1)

    Ms = 8.6e5
    A = 1.3e-11
    D = 1.23e-3

    function D_fun(i, j, k, dx, dy, dz)
        return D
    end

    sim = Sim(mesh)
    set_Ms(sim, Ms)
    init_m0(sim, m0_fun; norm=false)
    dmi = add_dmi(sim, D)
    dmi2 = add_dmi(sim, D_fun; name="dmi2")

    NuMag.effective_field(sim, sim.spin, 0.0)
    f1 = Array(dmi.field)
    f2 = Array(dmi2.field)

    @test isapprox(f1, f2, atol=1e-8)

    if isa(sim.spin, Array)
        NuMag.effective_field_debug(dmi, sim, sim.spin, 0.0)
        @test isapprox(f1, dmi.field, atol=1e-10)
    end

    dx = 2e-9

    expected_fx = -1.0 * (D / dx / (NuMag.mu_0 * Ms)) * NuMag.cross_x(1, 0, 0, 4, 5, 6)
    expected_fy = -1.0 * (D / dx / (NuMag.mu_0 * Ms)) * NuMag.cross_y(1, 0, 0, 4, 5, 6)
    expected_fz = -1.0 * (D / dx / (NuMag.mu_0 * Ms)) * NuMag.cross_z(1, 0, 0, 4, 5, 6)

    @test isapprox(f1[1], expected_fx)
    @test isapprox(f1[2], expected_fy)
    @test isapprox(f1[3], expected_fz)
end

function test_interfacial_dmi()
    function m0_fun(i, j, k, dx, dy, dz)
        if i == 1
            return (1, 2, 3)
        elseif i == 2
            return (4, 5, 6)
        end
        return (1, 0, 0)
    end

    mesh = FDMesh(; dx=2e-9, nx=2, ny=1, nz=1)

    Ms = 8.6e5
    A = 1.3e-11
    D = 1.23e-3

    sim = Sim(mesh)
    set_Ms(sim, Ms)
    init_m0(sim, m0_fun; norm=false)
    dmi = add_dmi(sim, D; type="interfacial")

    NuMag.effective_field(sim, sim.spin, 0.0)
    f1 = Array(dmi.field)

    dx = 2e-9

    expected_fx = 1.0 * (D / dx / (NuMag.mu_0 * Ms)) * NuMag.cross_x(0, -1, 0, 4, 5, 6)
    expected_fy = 1.0 * (D / dx / (NuMag.mu_0 * Ms)) * NuMag.cross_y(0, -1, 0, 4, 5, 6)
    expected_fz = 1.0 * (D / dx / (NuMag.mu_0 * Ms)) * NuMag.cross_z(0, -1, 0, 4, 5, 6)

    @test isapprox(f1[1], expected_fx)
    @test isapprox(f1[2], expected_fy)
    @test isapprox(f1[3], expected_fz)
end

@testset "Test DMI CPU" begin
    set_backend("cpu")
    test_bulk_dmi()
    test_interfacial_dmi()
end

@testset "Test DMI CUDA" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        test_bulk_dmi()
        test_interfacial_dmi()
    end
end
