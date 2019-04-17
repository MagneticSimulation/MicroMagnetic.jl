using JuMag
using Test

function test_FDMesh(mesh)
	@test mesh.volume == mesh.dx*mesh.dy*mesh.dz
	nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
	for k = 1:nz, j = 1:ny, i=1:nx
	  id = JuMag.index(i,j,k, nx, ny, nz)
	  @test mesh.ngbs[1, id] == JuMag._x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
	  @test mesh.ngbs[2, id] == JuMag._x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
	  @test mesh.ngbs[3, id] == JuMag._y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
	  @test mesh.ngbs[4, id] == JuMag._y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
	  @test mesh.ngbs[5, id] == JuMag._z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
	  @test mesh.ngbs[6, id] == JuMag._z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
    end
end

mesh1 = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=4, nz=2)
mesh2 = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=4, nz=2, pbc="xyz")

for mesh in (mesh1, mesh2)
    test_FDMesh(mesh)
end
