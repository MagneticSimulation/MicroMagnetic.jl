using MicroMagnetic
using Test
include("test_utils.jl")

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

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    f1 = Array(dmi.field)
    f2 = Array(dmi2.field)

    rtol =  MicroMagnetic.Float[] == Float64 ? 1e-8 : 1e-6 
    @test isapprox(f1, f2, rtol=rtol)

    if isa(sim.spin, Array{Float64})
        MicroMagnetic.effective_field_debug(dmi, sim, sim.spin, 0.0)
        @test isapprox(f1, dmi.field, atol=1e-10)
    end

    dx = 2e-9

    expected_fx = -1.0 * (D / dx / (MicroMagnetic.mu_0 * Ms)) * MicroMagnetic.cross_x(1, 0, 0, 4, 5, 6)
    expected_fy = -1.0 * (D / dx / (MicroMagnetic.mu_0 * Ms)) * MicroMagnetic.cross_y(1, 0, 0, 4, 5, 6)
    expected_fz = -1.0 * (D / dx / (MicroMagnetic.mu_0 * Ms)) * MicroMagnetic.cross_z(1, 0, 0, 4, 5, 6)

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

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    f1 = Array(dmi.field)

    dx = 2e-9

    expected_fx = 1.0 * (D / dx / (MicroMagnetic.mu_0 * Ms)) * MicroMagnetic.cross_x(0, -1, 0, 4, 5, 6)
    expected_fy = 1.0 * (D / dx / (MicroMagnetic.mu_0 * Ms)) * MicroMagnetic.cross_y(0, -1, 0, 4, 5, 6)
    expected_fz = 1.0 * (D / dx / (MicroMagnetic.mu_0 * Ms)) * MicroMagnetic.cross_z(0, -1, 0, 4, 5, 6)

    @test isapprox(f1[1], expected_fx)
    @test isapprox(f1[2], expected_fy)
    @test isapprox(f1[3], expected_fz)
end

@using_gpu()
test_functions("DMI", test_bulk_dmi, test_interfacial_dmi)
