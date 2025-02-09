using MicroMagnetic
using Test
using LinearAlgebra
using Enzyme

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function test_cubic_anis()
    mesh = CubicMesh(; nx=10, ny=1, nz=1)

    
    sim = Sim(mesh)

    mu_s = 2*mu_B
    set_mu_s(sim, mu_s)
    
    Kc = 1meV
    
    init_m0(sim, (0.6, 0.8, 0); norm=false)

    anis = add_cubic_anis(sim, Kc)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    @test isapprox(field[1], 1 / (mu_s) * 4 * Kc * 0.6^3)
    @test isapprox(field[2], 1 / (mu_s) * 4 * Kc * 0.8^3)
    @test isapprox(field[3], 1 / (mu_s) * 4 * Kc * 0^3)
    @test isapprox(energy[1], -Kc * (0.6^4 + 0.8^4))
end

function hexagonal_energy(m, K1, K2, K3)
    mx, my, mz = m[1], m[2], m[3]
    return K1*(1-mz*mz) + K2*(1-mz*mz)^2 + K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)
end

function test_hex_anis()
    mesh = CubicMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    mu_s = 2*mu_B
    set_mu_s(sim, mu_s)
    m0 = (0.7, -0.4, 1.2)
    init_m0(sim, m0; norm=false)

    K1, K2, K3 = 1.23e2*meV, 3.7e3*meV, 6.9e2*meV
    anis = add_hex_anis(sim, K1=K1, K2=K2, K3=K3)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    gd = Enzyme.gradient(Forward, hexagonal_energy, m0, Const(K1), Const(K2), Const(K3))
    expected = - collect(gd[1]) ./ mu_s

    @test isapprox(field[1:3], expected)
    @test isapprox(energy[1], hexagonal_energy(m0, K1, K2, K3), rtol=1e-5)
end

@using_gpu()
test_functions("Anisotropy", test_cubic_anis, test_hex_anis)