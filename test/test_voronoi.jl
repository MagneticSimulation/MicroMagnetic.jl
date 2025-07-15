using MicroMagnetic
using Test


if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_voronoi()
    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=2)

    grain_ids, gb_mask, points = voronoi(mesh, min_dist=30, seed=10000)
    
    #println(length(points))
    #println(grain_ids)
    @test size(grain_ids) == (100, 100)
    @test size(gb_mask) == (100, 100)
    @test all(grain_ids .> 0)  
    @test all(grain_ids .<= 35)  
    @test length(points) == 35
end

function test_plot_voronoi()
    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=2)

    grain_ids, gb_mask, points = voronoi(mesh, min_dist=30, seed=10000)

    plot_voronoi(grain_ids, points, dx=2, dy=2)
end


test_functions("Voronoi", test_voronoi)

#using CairoMakie
#test_plot_voronoi()