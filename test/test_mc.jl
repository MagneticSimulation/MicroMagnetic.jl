
using MicroMagnetic
using Test

#@using_gpu()

function test_monte_carlo(T=300)
    mesh = CubicMesh(; nx=1, ny=1, nz=1, pbc="x")
    sim = MonteCarlo(mesh; name="mc")
    init_m0(sim, (0, 0, 1))

    add_exch(sim; J=300 * k_B)
    e1 = MicroMagnetic.compute_system_energy(sim)
    @test isapprox(e1,-600)

    add_dmi(sim; D=0)
    add_zeeman(sim; Hx=0, Hy=0, Hz=100*k_B)
    e1 = MicroMagnetic.compute_system_energy(sim)
    @test isapprox(e1,-600-100)

    add_anis(sim; Ku=30*k_B, Kc=0)
    e1 = MicroMagnetic.compute_system_energy(sim)
    @test isapprox(e1,-600-100-30)
    
    sim.T = 100000
    run_sim(sim; maxsteps=1, save_vtk_every=-1, save_m_every=-1)

    MicroMagnetic.compute_system_energy(sim)
    
    t = MicroMagnetic.average_m(sim)
    return t
end

@testset "Test MC CPU" begin
    set_backend("cpu")
    test_monte_carlo()
end

@testset "Test MC CUDA" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        set_backend("cuda")
        test_monte_carlo()
    end
end
