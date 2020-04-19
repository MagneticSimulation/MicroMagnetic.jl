using JuMag
using Test

function test_geo(;gpu=false)
    Nz = 10
    mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=30, ny=30, nz=Nz)
    geo = JuMag.create_cylinder(mesh, JuMag.ex)
    @test isapprox(geo.xc,30e-9)
    @test isapprox(geo.yc,30e-9)
    @test isapprox(geo.zc,10e-9)
    geo = JuMag.create_cylinder(mesh, JuMag.ez)
    r = 30e-9

    ratio = sum(geo.shape)/Nz*mesh.dx*mesh.dy/(pi*r*r)
    @test ratio > 0.997
    println(ratio)
end

function test_set_Ms_cylindrical(;gpu=false)
    mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=5, ny=5, nz=1)
    Ms = 8.6e5
    sim = Sim(mesh)
    set_Ms_cylindrical(sim, Ms)
    @test sim.Ms[1] == 0
    @test sim.Ms[12] > 8e5
    println(sim.Ms[12])
end

@testset "Geometry" begin
    test_geo(gpu=false)
    test_set_Ms_cylindrical(gpu=false)
end

if JuMag._cuda_available.x
  JuMag.cuda_using_double()
  @testset "Geometry" begin
      test_geo(gpu=true)
      test_set_Ms_cylindrical(gpu=true)
  end
end
