using MicroMagnetic
using Test
using LinearAlgebra
using Enzyme

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_anis()
    mesh = FDMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)

    m = (1.1, 2.3, 4.3)
    init_m0(sim, m; norm=false)

    Ku = 1e5
    axis = (0, 0.6, 0.8)
    anis = add_anis(sim, Ku; axis=axis)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    @test isapprox(field[1], 0)
    @test isapprox(field[2], 2 * Ku / (MicroMagnetic.mu_0 * Ms) * dot(m, axis) * axis[2])
    @test isapprox(field[3], 2 * Ku / (MicroMagnetic.mu_0 * Ms) * dot(m, axis) * axis[3])
    @test isapprox(energy[10], Ku * (1 - dot(m, axis)^2) * 1e-27)
end

function test_cubic_anis()
    mesh = FDMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    Ms = 8.6e5
    Kc = 1e3
    set_Ms(sim, Ms)
    init_m0(sim, (0.6, 0.8, 0); norm=false)

    anis = add_cubic_anis(sim, Kc)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    @test isapprox(field[1], 1 / (MicroMagnetic.mu_0 * Ms) * 4 * Kc * 0.6^3)
    @test isapprox(field[2], 1 / (MicroMagnetic.mu_0 * Ms) * 4 * Kc * 0.8^3)
    @test isapprox(field[3], 1 / (MicroMagnetic.mu_0 * Ms) * 4 * Kc * 0^3)
    @test isapprox(energy[1], -Kc * (0.6^4 + 0.8^4) * 1e-27)
end

function hexagonal_energy(m, K1, K2, K3)
    mx, my, mz = m[1], m[2], m[3]
    return K1*(1-mz*mz) + K2*(1-mz*mz)^2 + K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)
end

function test_hex_anis()
    mesh = FDMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)
    m0 = (0.7, -0.4, 1.2)
    init_m0(sim, m0; norm=false)

    K1, K2, K3 = 1.23e2, 3.7e3, 6.9e2
    anis = add_hex_anis(sim, K1=K1, K2=K2, K3=K3)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    gd = Enzyme.gradient(Forward, hexagonal_energy, m0, Const(K1), Const(K2), Const(K3))
    expected = - collect(gd[1]) ./ (MicroMagnetic.mu_0*Ms)

    @test isapprox(field[1:3], expected)
    @test isapprox(energy[1]*1e27, hexagonal_energy(m0, K1, K2, K3), rtol=1e-5)
end

@using_gpu()
test_functions("Anisotropy", test_anis, test_cubic_anis, test_hex_anis)