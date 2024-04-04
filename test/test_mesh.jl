using JuMag
using Test

function test_FDMesh()
	mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=4, nz=2, pbc="xz")
	@test mesh.volume == mesh.dx*mesh.dy*mesh.dz
	nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
	ngbs = Array(mesh.ngbs)
	for k = [1, nz], j = [1, ny], i=1:nx
	  id = JuMag.index(i,j, k, nx, ny, nz)
	  @test ngbs[1, id] == JuMag._x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
	  @test ngbs[2, id] == JuMag._x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
	  @test ngbs[3, id] == JuMag._y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
	  @test ngbs[4, id] == JuMag._y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
	  @test ngbs[5, id] == JuMag._z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
	  @test ngbs[6, id] == JuMag._z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
    end
end

