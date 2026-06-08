using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_tensor()
    println("\n=== Test: Tensor Model ===")
    
    Ms = 8.6e5
    lambda_s = 30e-6
    mu0 = 4π * 1e-7
    
    mesh = FDMesh(; nx=10, ny=1, nz=1)
    sim = Sim(mesh)
    set_Ms(sim, Ms)
    
    # Create stress tensor: σxx = σ, others = 0
    sigma = zeros(10, 1, 1, 6)
    sigma[:, :, :, 1] .= 1e8  # σxx = 100 MPa
    println("  DEBUG: input sigma[1,1,1,:] = ", sigma[1,1,1,:])
    
    # Test: m = (1, 0, 0)
    init_m0(sim, (1, 0, 0))
    me = add_mel(sim; model=:tensor, lambda_s=lambda_s, stress_or_strain=sigma)
    
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    
    println("  DEBUG: me.field_data = ", me.field_data[1:6])
    println("  DEBUG: mesh.volume = ", mesh.volume)
    
    # Energy: E = -1.5 * λs * σxx * mx² * volume * n_cells
    n_cells = mesh.nx * mesh.ny * mesh.nz
    expected_E = -1.5 * lambda_s * 1e8 * mesh.volume * n_cells
    @test isapprox(sum(me.energy), expected_E, atol=abs(expected_E) * 1e-5)
    
    # Field: Hx = 3*λs/(μ₀*Ms) * σxx * mx
    factor = 3 * lambda_s / (mu0 * Ms)
    expected_Hx = factor * 1e8 * 1.0
    field = Array(me.field)
    @test isapprox(field[1], expected_Hx, atol=abs(expected_Hx) * 1e-5)
    
    println("  ✓ Energy and field correct for σxx tensor")
    
    return true
end

function test_cubic()
    println("\n=== Test: Cubic Model ===")
    
    Ms = 8.6e5
    B1 = 1e6
    B2 = 0.5e6
    mu0 = 4π * 1e-7
    
    mesh = FDMesh(; nx=10, ny=1, nz=1)
    sim = Sim(mesh)
    set_Ms(sim, Ms)
    
    # Create strain: εxx = ε, εyy = ε, others = 0
    eps = 1e-3
    strain = zeros(10, 1, 1, 6)
    strain[:, :, :, 1] .= eps  # εxx
    strain[:, :, :, 2] .= eps  # εyy
    
    # Test: m = (1, 0, 0)
    init_m0(sim, (1, 0, 0))
    me = add_mel(sim; model=:cubic, B1=B1, B2=B2, stress_or_strain=strain)
    
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    
    # Energy: E = B1 * εxx * mx² * volume * n_cells
    n_cells = 10
    expected_E = B1 * eps * mesh.volume * n_cells
    @test isapprox(sum(me.energy), expected_E, atol=abs(expected_E) * 1e-5)
    
    # Field: Hx = -2/(μ₀*Ms) * B1 * εxx * mx
    factor = 2 / (mu0 * Ms)
    expected_Hx = -factor * B1 * eps
    field = Array(me.field)
    @test isapprox(field[1], expected_Hx, atol=abs(expected_Hx) * 1e-5)
    @test isapprox(field[2], 0.0, atol=1e-10)
    @test isapprox(field[3], 0.0, atol=1e-10)
    
    println("  ✓ Energy and field correct for εxx strain")
    
    # Test: shear strain with B2
    strain2 = zeros(10, 1, 1, 6)
    strain2[:, :, :, 4] .= eps  # εxy
    
    sim2 = Sim(mesh)
    set_Ms(sim2, Ms)
    init_m0(sim2, (1/sqrt(2), 1/sqrt(2), 0))
    me2 = add_mel(sim2; model=:cubic, B1=0, B2=B2, stress_or_strain=strain2)
    
    MicroMagnetic.effective_field(sim2, sim2.spin, 0.0)
    
    # Energy: E = 2*B2 * εxy * mx*my * volume * n_cells = B2 * ε * volume * n_cells
    expected_E2 = B2 * eps * mesh.volume * n_cells
    @test isapprox(sum(me2.energy), expected_E2, atol=abs(expected_E2) * 1e-5)
    
    # Field: Hx = -2*B2/(μ₀*Ms) * εxy * my, Hy = -2*B2/(μ₀*Ms) * εxy * mx
    factor2 = 2 / (mu0 * Ms)
    expected_H = -factor2 * B2 * eps / sqrt(2)
    field2 = Array(me2.field)
    @test isapprox(field2[1], expected_H, atol=abs(expected_H) * 1e-5)
    @test isapprox(field2[2], expected_H, atol=abs(expected_H) * 1e-5)
    
    println("  ✓ B2 term with factor 2 correct")
    
    return true
end

@using_gpu()
test_functions("Magnetoelastic", test_tensor, test_cubic)
