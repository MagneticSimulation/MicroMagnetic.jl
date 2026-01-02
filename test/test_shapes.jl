using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_shapes()

    p1 = Plane(point=(10,0,0), normal=(1, 0, 0))
    @test MicroMagnetic.inside(p1, (9.9,0,0)) == false
    @test MicroMagnetic.inside(p1, (10.1,0,0)) == true
    @test MicroMagnetic.inside(p1, (Inf,0,0)) == true

    c1 = Cylinder(radius=50e-9, normal=(0,0,1))
    @test MicroMagnetic.inside(c1, (20e-9,30e-9,100)) == true
    @test MicroMagnetic.inside(c1, (50.1e-9,0,0)) == false

    s1 = Sphere(radius = 30e-9, center=(30e-9, 30e-9, 0))
    @test MicroMagnetic.inside(s1, (20e-9,30e-9,100)) == false
    @test MicroMagnetic.inside(s1, (20e-9,30e-9,0)) == true

    b1 = Box(sides = (50e-9, 50e-9, Inf), theta=pi/4)
    @test MicroMagnetic.inside(b1, (50e-9*sqrt(2)/2*0.9999,0,100)) == true
    @test MicroMagnetic.inside(b1, (50e-9,0,0)) == false

    t1 = Torus(R = 60e-9, r=20e-9)
    @test MicroMagnetic.inside(t1, (50e-9,0,0)) == true
    @test MicroMagnetic.inside(t1, (0,0,0)) == false

    c2 = t1 + p1 * s1
    @test MicroMagnetic.inside(c2, (50e-9,0,0)) == true
    @test MicroMagnetic.inside(c2, (0,0,0)) == false

    c3 = t1 - (p1 * s1)
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=50)
    save_vtk(mesh, c3, "shape")
    
end

@using_gpu()
test_functions("Geometry", test_shapes)