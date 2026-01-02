using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_sim()
    #Test mesh
    mesh = FDMesh(; dx=1.1e-9, nx=10, ny=2)

    sim = Sim(mesh; driver="LLG")
    set_Ms(sim, 1.0)

    @test sim.n_total == 20

    init_m0(sim, (1.0, 1.0, 0))
    #println(sim.spin)
    spin = Array(sim.spin)

    @test isapprox(spin[1], 0.5 * 2^0.5)
    @test isapprox(spin[2], 0.5 * 2^0.5)


    m0 = Float32[i for i = 1:60]
    init_m0(sim, m0, norm=false)
    @test isapprox(Array(sim.spin), m0)

    set_Ms(sim, 8.6e5)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.21e5
    sim.driver.precession = false

    add_exch(sim, 1.3e-11)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    @info("test_sim() passed!")
    #println(sim.spin)
end

function test_region_map()
    # Basic region_map functionality
    Ms_map = region_map(1 => 8.6e5)
    @test Ms_map(1) == 8.6e5
    @test Ms_map(2) == 0.0
    @test Ms_map(0) == 0.0
    
    # Multiple regions
    exch_map = region_map(
        1 => 1.3e-11,
        2 => 0.8e-11,
        3 => 0.5e-11
    )
    @test exch_map(1) == 1.3e-11
    @test exch_map(2) == 0.8e-11
    @test exch_map(3) == 0.5e-11
    @test exch_map(4) == 0.0
    
    # With custom default value
    anis_map = region_map(
        -1 => 5e5,
        2 => 2e5,
        default=1e5
    )
    @test anis_map(-1) == 5e5
    @test anis_map(2) == 2e5
    @test anis_map(3) == 1e5
    @test anis_map(0) == 1e5
    
    # Integration with set_region and set_Ms (basic)
    mesh = FDMesh(nx=10, ny=10, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
    circle = Cylinder(center=(0, 0, 0), radius=20e-9, height=10e-9, normal=(0, 0, 1))
    set_region(mesh, circle, 1)
    
    # Verify set_region worked correctly
    regions = Array(mesh.regions)
    @test any(region_id -> region_id == 1, regions)
    @test any(region_id -> region_id == 0, regions)
    
    # Test region_map integration
    sim = Sim(mesh; driver="LLG")
    set_Ms(sim, region_map(1 => 8.6e5, default=0.0))
    Ms_array = Array(sim.mu0_Ms)
    @test any(ms -> ms > 0.0, Ms_array)  # Non-zero Ms inside region 1
    @test any(ms -> ms ≈ 0.0, Ms_array)  # Zero Ms outside region 1
end

@using_gpu()
test_functions("Sim", test_sim, test_region_map)





