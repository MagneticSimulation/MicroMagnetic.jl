using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_df_torque(integrator="DormandPrince")
    # Test mesh
    mesh = FDMesh(; nx=1, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.tol = 1e-8

    # Add Zeeman field for reference
    add_zeeman(sim, (0, 0, 1e5))

    # Test parameters
    aj = 1e9  
    bj = 5e8 
    p = (0.1, 0.2, 0.3)  # polarization direction

    # Add DFT torque
    t = add_sot(sim, aj, bj, p; name="dft_torque")

    # Initialize magnetization
    init_m0(sim, (-0.3, 0.4, 0.8), norm=false)
    
    # Calculate effective field
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    
    # Expected effective field calculation based on the formula:
    # H_stt = (1/γ) * (a_J * (m × p) + b_J * p)
    # where m = (-0.3, 0.4, 0.8), p = (0.1, 0.2, 0.3), γ = 2.21e5
    
    # Calculate m × p
    m_cross_p = MicroMagnetic.cross([-0.3, 0.4, 0.8], [0.1, 0.2, 0.3])
    torque_field = aj .* m_cross_p .+ bj .* [0.1, 0.2, 0.3]
    expected = torque_field ./ 2.21e5

    # Test if the calculated field matches expected
    @test isapprox(Array(t.field), expected, rtol=1e-6)
    
    return true
end

function test_df_torque_atomistic(integrator="DormandPrince")
    # Test mesh for atomistic simulation
    mesh = CubicMesh(; nx=1, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_mu_s(sim, 2*mu_B)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5/mu_0  # Note: atomistic gamma is divided by mu_0
    sim.driver.integrator.tol = 1e-8

    # Add Zeeman field for reference
    add_zeeman(sim, (0, 0, 1e5))

    # Test parameters (same as continuum)
    aj = 1e9
    bj = 5e8
    p = (0.1, 0.2, 0.3)  # polarization direction

    t = add_sot(sim, aj, bj, p; name="dft_torque_atomistic")

    # Initialize magnetization
    init_m0(sim, (-0.3, 0.4, 0.8), norm=false)
    
    # Calculate effective field
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    
    # Expected calculation (same physics, different gamma)
    m_cross_p = MicroMagnetic.cross([-0.3, 0.4, 0.8], [0.1, 0.2, 0.3])
    torque_field = aj .* m_cross_p .+ bj .* [0.1, 0.2, 0.3]
    expected = torque_field ./ (2.21e5/mu_0)  # Note: atomistic gamma is divided by mu_0
    
    # Test if the calculated field matches expected
    @test isapprox(Array(t.field), expected, rtol=1e-6)
    
    return true
end

function test_df_torque_functional_aj(integrator="DormandPrince")
    # Test with functional aj parameter
    mesh = FDMesh(; nx=2, ny=2, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5

    # Functional aj that depends on position
    function aj_func(i, j, k, dx, dy, dz)
        x = i * dx
        return 1e9 * (1.0 + 0.1 * x / 1e-9)  # Linear gradient in x-direction
    end

    bj = 5e8
    p = (0.1, 0.2, 0.3)

    t = add_sot(sim, aj_func, bj, p; name="dft_torque_func")

    init_m0(sim, (-0.3, 0.4, 0.8), norm=false)
    
    # Calculate effective field
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    
    # For functional aj, we can't easily compute exact expected values
    # but we can test that the field is non-zero and has expected dimensions
    @test size(t.field) == (12, )
    @test any(t.field .!= 0)  # At least some components should be non-zero
    
    return true
end

function test_df_torque_zero_polarization(integrator="DormandPrince")
    # Test edge case: zero polarization vector
    mesh = FDMesh(; nx=1, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.gamma = 2.21e5

    aj = 1e9
    bj = 5e8
    p = (0.0, 0.0, 0.0)  # Zero polarization

    # This should work without errors
    t = add_sot(sim, aj, bj, p; name="dft_torque_zero_p")

    init_m0(sim, (-0.3, 0.4, 0.8), norm=false)
    
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    
    # With zero polarization, the torque field should be zero
    @test all(isapprox.(Array(t.field), 0.0, atol=1e-12))
    
    return true
end

@using_gpu()
test_functions("DFTorque", 
    test_df_torque, 
    test_df_torque_atomistic, 
    test_df_torque_functional_aj,
    test_df_torque_zero_polarization
)
