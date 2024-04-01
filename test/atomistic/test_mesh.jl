using JuMag
using Test

function test_CylindricalTubeMeshGPU(mesh)
	@test mesh.n_total == mesh.nr*mesh.nz
	nr, nz = mesh.nr, mesh.nz
	ngbs = Array(mesh.ngbs)
	for k = 1:nz, i=1:mesh.nr
	  id = JuMag.index(i, 1, k, mesh.nr, 1, nz)
	  @test ngbs[1, id] == JuMag._x_minus_one(i, id, mesh.nr, 1, nz, true)
	  @test ngbs[2, id] == JuMag._x_plus_one(i, id, mesh.nr, 1, nz, true)
	  @test ngbs[3, id] == JuMag._z_minus_one(k, id, mesh.nr, 1, nz, mesh.zperiodic)
	  @test ngbs[4, id] == JuMag._z_plus_one(k, id, mesh.nr, 1, nz, mesh.zperiodic)
    end
end


mesh1  = CylindricalTubeMeshGPU(nz=2, nr=2, R=3, dz=2)
mesh2  = CylindricalTubeMeshGPU(nz=3, nr=4, R=3, dz=2, pbc="z")

for mesh in (mesh1, mesh2)
    test_CylindricalTubeMeshGPU(mesh)
end
