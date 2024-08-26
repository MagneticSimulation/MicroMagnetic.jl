using MicroMagnetic
using Test

@testset "Test ovfs cuda" begin
    m = [0.6, 0.8, 0, 0.6, 0.8, 0]
    mesh = FDMesh(; nx=2, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)
    MicroMagnetic.cuda_using_double(true)
    mesh = FDMeshGPU(; nx=2, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)
    sim = Sim(mesh; name="double_saving")
    set_Ms(sim, 8.0e5)
    init_m0(sim, (0.6, 0.8, 0))
    save_ovf(sim, "double_saved_double"; type=Float64)
    save_ovf(sim, "double_saved_single"; type=Float32)
    save_ovf(sim, "double_saved_string"; type=String)

    MicroMagnetic.cuda_using_double(false)
    sim = Sim(mesh; name="single_saving")
    set_Ms(sim, 8.0e5)
    init_m0(sim, (0.6, 0.8, 0))
    save_ovf(sim, "single_saved_double"; type=Float64)
    save_ovf(sim, "single_saved_single"; type=Float32)
    save_ovf(sim, "single_saved_string"; type=String)

    testing_files = ["double_saved_double", "double_saved_single", "double_saved_string",
                     "single_saved_double", "single_saved_single", "single_saved_string"]

    function run_test()
        sim = Sim(mesh)
        set_Ms(sim, 8.0e5)
        for f in testing_files
            read_ovf(sim, f)
            spin = Array(sim.spin)
            @test isapprox(spin, m, atol=1e-6)
            spin = Array(sim.prespin)
            @test isapprox(spin, m, atol=1e-6)
            relax(sim; max_steps=20, save_m_every=-1)
        end
    end

    MicroMagnetic.cuda_using_double(true)
    run_test()

    MicroMagnetic.cuda_using_double(false)
    run_test()

    for f in testing_files
        rm(f * ".ovf")
    end
end
