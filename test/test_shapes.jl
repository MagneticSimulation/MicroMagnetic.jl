using JuMag
using Test

function test_shapes()

    p1 = Plane(point=(10,0,0), normal=(1, 0, 0))
    @test JuMag.point_inside_shape((9.9,0,0), p1) == false
    @test JuMag.point_inside_shape((10.1,0,0), p1) == true
    @test JuMag.point_inside_shape((Inf,0,0), p1) == true

    c1 = Cylinder(radius=50e-9, normal=(0,0,1))
    @test JuMag.point_inside_shape((20e-9,30e-9,100), c1) == true
    @test JuMag.point_inside_shape((50.1e-9,0,0), c1) == false

    s1 = Sphere(radius = 30e-9, center=(30e-9, 30e-9, 0))
    @test JuMag.point_inside_shape((20e-9,30e-9,100), s1) == false
    @test JuMag.point_inside_shape((20e-9,30e-9,0), s1) == true

    b1 = Box(sides = (50e-9, 50e-9, Inf), theta=pi/4)
    @test JuMag.point_inside_shape((50e-9*sqrt(2)/2*0.9999,0,100), b1) == true
    @test JuMag.point_inside_shape((50e-9,0,0), b1) == false

    t1 = Torus(R = 60e-9, r=20e-9)
    @test JuMag.point_inside_shape((50e-9,0,0), t1) == true
    @test JuMag.point_inside_shape((0,0,0), t1) == false

    c2 = t1 + p1 * s1
    @test JuMag.point_inside_shape((50e-9,0,0), c2) == true
    @test JuMag.point_inside_shape((0,0,0), c2) == false

    c3 = t1 - (p1 * s1)
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=50)
    save_vtk(mesh, c3, "shape")
    
end

@testset "Geometry" begin
    test_shapes()
end

