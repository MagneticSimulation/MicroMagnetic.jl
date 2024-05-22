using MicroMagnetic
using Test

omega = Array([0.1,0.2,0.3,0,0,0,0.2,0.3,0.4])
spin = Array([0.6,0.8,0,0.6,0.8,0,0.6,0,0.8])
spin_next = zeros(Float64,9)

MicroMagnetic.omega_to_spin(omega, spin, spin_next, 3)

expected = Array([0.33816425120772947, 0.9410628019323674, -0.006763285024154589,
                 0.6, 0.8, 0.0, 0.7836829836829837, 0.1361305361305361, 0.6060606060606062])
@test isapprox(spin_next, expected)


mesh1 = FDMesh(nx=3, ny=5, nz=7)
function test_init_scalar(mesh)
	function index_wrapper(i,j,k,dx,dy,dz)
		return MicroMagnetic.index(i,j,k,mesh.nx, mesh.ny, mesh.nz)
	end
	n_total = mesh.n_total
	v = zeros(Float64, n_total)
	MicroMagnetic.init_scalar!(v, mesh, index_wrapper)
	v3d = reshape(v, mesh.nx, mesh.ny, mesh.nz)
	id = 0.0
    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        id += 1.0
        @test v3d[i,j,k] == id
    end
end

test_init_scalar(mesh1)
