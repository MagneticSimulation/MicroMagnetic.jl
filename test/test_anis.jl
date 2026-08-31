using MicroMagnetic
using Test
using LinearAlgebra

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


function test_spatial_anis()
    mesh = FDMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)

    m = (1.1, 2.3, 4.3)
    init_m0(sim, m; norm=false)

    Ku = 1e5
    axis_fun = (i,j,k,dx,dy,dz) -> (0, 0.6, 0.8)
    axis = (0, 0.6, 0.8)
    anis = add_anis(sim, Ku; axis=axis_fun)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    @test isapprox(field[1], 0)
    @test isapprox(field[2], 2 * Ku / (MicroMagnetic.mu_0 * Ms) * dot(m, axis) * axis[2])
    @test isapprox(field[3], 2 * Ku / (MicroMagnetic.mu_0 * Ms) * dot(m, axis) * axis[3])
    @test isapprox(energy[10], Ku * (1 - dot(m, axis)^2) * 1e-27)
    # a second site locks the per-spin component layout of the spatial axis
    @test isapprox(field[4], field[1])
    @test isapprox(field[5], field[2])
    @test isapprox(field[6], field[3])
    @test isapprox(energy[2], energy[10])
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

# Analytic gradient of hexagonal_energy (replaces Enzyme.gradient)
#   H = - (dE/dm) / (mu_0 Ms)  for FD  /  - (dE/dm) / mu_s for atomistic
function hexagonal_energy_gradient(m, K1, K2, K3)
    mx, my, mz = m[1], m[2], m[3]
    mx2, my2, mz2 = mx*mx, my*my, mz*mz
    dEmx = K3 * (6 * mx^5 - 60 * mx^3 * my2 + 30 * mx * my2^2)
    dEmy = K3 * (-30 * mx2^2 * my + 60 * mx2 * my^3 - 6 * my^5)
    dEmz = -2 * K1 * mz - 4 * K2 * (1 - mz2) * mz
    return (dEmx, dEmy, dEmz)
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

    dEmx, dEmy, dEmz = hexagonal_energy_gradient(m0, K1, K2, K3)
    expected = [ -dEmx, -dEmy, -dEmz ] ./ (MicroMagnetic.mu_0 * Ms)

    @test isapprox(field[1:3], expected)
    @test isapprox(energy[1]*1e27, hexagonal_energy(m0, K1, K2, K3), rtol=1e-5)
end

@using_gpu()
test_functions("Anisotropy", test_anis, test_spatial_anis, test_cubic_anis, test_hex_anis; platforms=["CPU"])
