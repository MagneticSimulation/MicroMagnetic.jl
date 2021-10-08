using JuMag
using Test

function test_geo(;gpu=false)
    Nz = 10
    r = 30e-9
    create_mesh = gpu ? FDMeshGPU : FDMesh
    mesh =  create_mesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=30, ny=30, nz=Nz)
    geo = JuMag.create_cylinder(mesh, JuMag.ex,r1=r,r2=r)

    @test isapprox(geo.xc,30e-9)
    @test isapprox(geo.yc,30e-9)
    @test isapprox(geo.zc,10e-9)
    ratio = sum(geo.shape)/Nz*mesh.dx*mesh.dy/(pi*r*r)
    @test ratio > 0.997
    println(ratio)
end

function test_set_Ms_cylindrical(;gpu=false)
    create_mesh = gpu ? FDMeshGPU : FDMesh
    mesh =  create_mesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=5, ny=5, nz=1)
    Ms = 8.6e5
    sim = Sim(mesh)
    set_Ms_cylindrical(sim, Ms, r1=5e-9, r2=5e-9)
    ms = Array(sim.Ms)
    @test ms[1] == 0
    @test ms[12] > 8e5
end

function test_functions_in_geometry(;gpu=false)
    create_mesh = gpu ? FDMeshGPU : FDMesh
    mesh =  create_mesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=5, ny=5, nz=5)
    sim = Sim(mesh)
    box = create_box(mesh,x1=0,y1=0,z1=0,x2=4e-9,y2=10e-9,z2=10e-9)
    set_Ms(sim,box,1e5)
    exch = add_exch(sim, box, 1e-12)
    anis = add_anis(sim, box, 1e5)
    ms = reshape(Array(sim.Ms), (5,5,5))
    As = reshape(Array(exch.A), (5,5,5))
    Ku = reshape(Array(anis.Ku), (5,5,5))
    for i = 1:5, j =1, k=1
        if i<=2
            @test ms[i,j,k] == 1e5
            @test As[i,j,k] == 1e-12
            @test Ku[i,j,k] == 1e5
        else
            @test ms[i,j,k] == 0
            @test As[i,j,k] == 0
            @test Ku[i,j,k] == 0
        end
    end
end

@testset "Geometry" begin
    test_geo(gpu=false)
    test_set_Ms_cylindrical(gpu=false)
    test_functions_in_geometry(;gpu=false)
end

if JuMag._cuda_available.x
  JuMag.cuda_using_double()
  @testset "Geometry" begin
      test_geo(gpu=true)
      test_set_Ms_cylindrical(gpu=true)
      test_functions_in_geometry(;gpu=true)
  end
end
