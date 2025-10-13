using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function spatial_alpha(i, j, k, dx, dy, dz)
    if i < 10
        return 0.01
    else
        return 0.99
    end
end

function test_spatial_llg()
    mesh = FDMesh(; nx=100, ny=1, nz=1, dx=2e-9, dy=2e-9, dz=1e-9)
    sim = Sim(mesh; name="test_alpha", driver="SpatialLLG")
    
    set_Ms(sim, 8.6e5)

    set_alpha(sim, (i, j, k, dx, dy, dz) -> 0.1)
    alpha = Array(sim.driver.alpha)
    @test maximum(alpha) == minimum(alpha)

    alpha_array = rand(sim.n_total) * 0.2 .+ 0.01
    set_alpha(sim, alpha_array)

    @test maximum(alpha) <= 0.21
    @test minimum(alpha) >= 0.01

    set_alpha(sim, spatial_alpha)
    alpha = Array(sim.driver.alpha)

    @test alpha[9] == 0.01
    @test alpha[10] == 0.99
end

@using_gpu()
test_functions("SpatialLLG", test_spatial_llg, precisions=[Float64])
