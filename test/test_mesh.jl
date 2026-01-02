using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_FDMesh()
    T = MicroMagnetic.Float[]
    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=4, nz=2, pbc="xz")
    @test isapprox(mesh.volume, 8e-27, rtol=1e-7)
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    ngbs = Array(mesh.ngbs)
    for k in [1, nz], j in [1, ny], i in 1:nx
        id = MicroMagnetic.index(i, j, k, nx, ny, nz)
        @test ngbs[1, id] == MicroMagnetic._x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        @test ngbs[2, id] == MicroMagnetic._x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        @test ngbs[3, id] == MicroMagnetic._y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        @test ngbs[4, id] == MicroMagnetic._y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        @test ngbs[5, id] == MicroMagnetic._z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        @test ngbs[6, id] == MicroMagnetic._z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
    end

    Nx = 50
    mesh = FDMesh(; dx=2e-9, nx=Nx, ny=1, nz=1, pbc="x")
    @test mesh.dx == T(2e-9)
    @test mesh.nx == Nx
    @test mesh.x0 + 50e-9 < eps()
    ngbs = Array(mesh.ngbs)
    @test ngbs[1, 1] == Nx
    @test ngbs[1, Nx] == Nx - 1
    @test ngbs[2, 1] == 2
    @test ngbs[2, Nx] == 1

    mesh = FDMesh(; dx=1.1e-9, nx=10, ny=2)
    @test isapprox(mesh.dx * 1.0, 1.1e-9, rtol=1e-7)
    @test mesh.nx == 10
    @test isapprox(mesh.volume * 1.0, 1.1e-27, rtol=1e-7)
end

function test_set_region()
    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=10, ny=10, nz=1)

    sphere = Sphere(radius=10e-9, center=(10e-9, 10e-9, 0))
    
    set_region(mesh, sphere, 1)
    
    regions = Array(mesh.regions)
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = MicroMagnetic.index(i, j, k, nx, ny, nz)
        x = mesh.x0 + (i - 0.5)*dx
        y = mesh.y0 + (j - 0.5)*dy
        z = mesh.z0 + (k - 0.5)*dz
        
        if MicroMagnetic.inside(sphere, (x, y, z))
            @test regions[id] == 1
        else
            @test regions[id] == 0
        end
    end
    
end

@using_gpu()
test_functions("Mesh", test_FDMesh, test_set_region)