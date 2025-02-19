using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function test_CylindricalTubeMesh(mesh)
    @test mesh.n_total == mesh.nr * mesh.nz
    nr, nz = mesh.nr, mesh.nz
    ngbs = Array(mesh.ngbs)
    for k in 1:nz, i in 1:(mesh.nr)
        id = MicroMagnetic.index(i, 1, k, mesh.nr, 1, nz)
        @test ngbs[1, id] == MicroMagnetic._x_minus_one(i, id, mesh.nr, 1, nz, true)
        @test ngbs[2, id] == MicroMagnetic._x_plus_one(i, id, mesh.nr, 1, nz, true)
        @test ngbs[3, id] == MicroMagnetic._z_minus_one(k, id, mesh.nr, 1, nz, mesh.zperiodic)
        @test ngbs[4, id] == MicroMagnetic._z_plus_one(k, id, mesh.nr, 1, nz, mesh.zperiodic)
    end
end

function test_TriangularMesh()
    mesh = TriangularMesh(dx=1e-9, dz=1e-9, nx=3, ny=3, nz=1, pbc="xy")
    @test mesh.nz == 1
    @test mesh.n_ngbs == 8
    @test mesh.n_ngbs2 == 6
    @test size(mesh.ngbs) == (8, mesh.n_total)
    @test size(mesh.ngbs2) == (6, mesh.n_total)
end


function test_CubicMesh()
    mesh = CubicMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=2, nz=3, pbc="x")
    @test mesh.nx == 10
    @test mesh.n_ngbs == 6
    @test mesh.n_ngbs2 == 12
    @test mesh.n_ngbs3 == 8
    @test mesh.n_ngbs4 == 6

    MicroMagnetic.compute_2nd_ngbs(mesh)
    @test size(mesh.ngbs2) == (12, mesh.n_total)

    MicroMagnetic.compute_3rd_ngbs(mesh)
    @test size(mesh.ngbs3) == (8, mesh.n_total)

    MicroMagnetic.compute_4th_ngbs(mesh)
    @test size(mesh.ngbs4) == (6, mesh.n_total)
end


function test_all_meshes()
    mesh1 = CylindricalTubeMesh(; nz=2, nr=2, R=3, dz=2)
	mesh2 = CylindricalTubeMesh(; nz=3, nr=4, R=3, dz=2, pbc="z")
	for mesh in (mesh1, mesh2)
		test_CylindricalTubeMesh(mesh)
	end

    test_TriangularMesh()
    test_CubicMesh()
end

@using_gpu()
test_functions("AtomisticMeshes", test_all_meshes)