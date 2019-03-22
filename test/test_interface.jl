using SpinDynamics
using Test

mesh = FDMesh(dx=1.0, dy=1.0, dz=1.0, nx=100, ny=100, nz=1, unit_length=1e-9)

function circular_Ms(i,j,k,dx,dy,dz)
	if (i-50.5)^2 + (j-50.5)^2 <= 50^2
		return 8.6e5
	end
	return 0.0
end

set_Ms(mesh, circular_Ms)

@test mesh.Ms[1] == 0.0
@test mesh.Ms[10] == 0.0
@test mesh.Ms[50] == 8.6e5

@test sim.nxyz == 10000
