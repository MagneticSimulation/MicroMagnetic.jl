using MicroMagnetic
using Test
using LinearAlgebra
using Random

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function init_fun3(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end

function init_fun6(i, j, k, dx, dy, dz)
    x = i - 50.5
    y = j - 50.5
    r = (x^2 + y^2)^0.5
    if r < 20
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end


function test_init_m0()
    mesh = FDMesh(; nx=100, ny=100, nz=2)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)

    init_m0(sim, init_fun3)

    m1 = Array(sim.spin)

    init_m0(sim, init_fun6)
    m2 = Array(sim.spin)

    @test isapprox(m1, m2)

end

function test_init_m0_random()
    mesh = FDMesh(; nx=10, ny=10, nz=1)
    
    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)

    init_m0_random(sim, seed=42)
    spin1 = Array(sim.spin)
    
    init_m0_random(sim, seed=42)
    spin2 = Array(sim.spin)
    
    @test isapprox(spin1, spin2)
    
    @test length(unique(spin1)) > 3
    
    @test all(-1.0 .<= spin1 .<= 1.0)
end

function test_vortex()
    mesh = FDMesh(; nx=50, ny=50, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)
    
    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)

    # Test vortex function with init_m0
    init_m0(sim, vortex(center=(25e-9, 25e-9), R=10e-9, p=1, c=1))
    
    b = reshape(Array(sim.spin), (3, 50, 50, 1))
    mz = b[3, :, :, :]
    
    @test minimum(mz) ≈ 0 
    @test maximum(mz) ≈ 1 
end

function test_skyrmion()
    mesh = FDMesh(; nx=50, ny=50, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)
    
    sim = Sim(mesh)
    set_Ms(sim, 8e5)

    # Test skyrmion function with init_m0 (Néel skyrmion)
    init_m0(sim, skyrmion((4.5e-9, 4.5e-9), 10e-9, 1, 1, 0.0))
    
    b = reshape(Array(sim.spin), (3, 50, 50, 1))
    
    # Check center region (should be up along z-axis)
    center_spin = b[:, 30, 30, 1]
    @test isapprox(center_spin[3], 1.0; atol=0.1)
    
    # Check edge region (should be down along z-axis)
    edge_spin = b[:, 1, 1, 1]
    @test isapprox(edge_spin[3], -1.0; atol=0.1)
    
    # Test skyrmion lattice (Bloch type)
    init_m0(sim, skyrmion_lattice(25e-9, p=-1, c=1, type=:bloch))
    b_lattice_bloch = reshape(Array(sim.spin), (3, 50, 50, 1))
    
    # Test skyrmion lattice (Néel type)
    init_m0(sim, skyrmion_lattice(25e-9, p=-1, c=1, type=:neel))
    b_lattice_neel = reshape(Array(sim.spin), (3, 50, 50, 1))
    
    # Check that Bloch and Néel lattices have different magnetization patterns
    @test !isapprox(b_lattice_bloch, b_lattice_neel; atol=0.1)
    
end

@using_gpu()
test_functions("init_m0", test_init_m0, test_init_m0_random, test_vortex, test_skyrmion)