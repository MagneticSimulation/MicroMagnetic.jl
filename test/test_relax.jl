using MicroMag
using Test

function init_dw(i, j, k, dx, dy, dz)
    if i < 240
        return (1, 0.1, 0)
    elseif i < 260
        return (0, 1, 0.1)
    else
        return (-1, 0.1, 0)
    end
end

function relax_system(; driver="LLG")
    mesh = FDMesh(; nx=500, ny=1, nz=11, dx=2e-9, dy=2e-9, dz=1e-9)
    sim = Sim(mesh; name="relax_llg", driver=driver)
    if driver == "LLG"
        sim.driver.precession = false
        sim.driver.alpha = 0.5
    end
    set_Ms(sim, 8.6e5)

    add_exch(sim, 1.3e-11)
    add_anis(sim, 1e5; axis=(1, 0, 0))

    init_m0(sim, init_dw)
    relax(sim; maxsteps=2000, stopping_dmdt=0.1)
    m = MicroMag.average_m(sim)
    @test abs(m[1]) < 0.01
    @test abs(m[2]) < 0.04
    @test abs(m[3]) < 0.01
end

@testset "Test Relax CPU" begin
    set_backend("cpu")
    relax_system(; driver="LLG")
    relax_system(; driver="SD")
end

@testset "Test Relax CUDA" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        relax_system(; driver="LLG")
        relax_system(; driver="SD")
    end
end
