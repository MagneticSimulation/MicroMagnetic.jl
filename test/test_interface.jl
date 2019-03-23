using SpinDynamics
using Test

mesh = FDMesh(dx=1.0, dy=1.0, dz=1.0, nx=100, ny=100, nz=1, unit_length=1e-9, pbc="xy")

function circular_Ms(i,j,k,dx,dy,dz)
	if (i-50.5)^2 + (j-50.5)^2 <= 50^2
		return 8.6e5
	end
	return 0.0
end

function m0_fun(i,j,k,dx,dy,dz)
  L = 50*dx
  x = i*dx
  return sin(2*pi*x/L), sin(2*pi*x/L+1.2), sin(2*pi*x/L+2.3)
end

function test_energy_minimization(mesh)
	sim = Sim(mesh, driver="SDM", name="sim")
	@test set_Ms(sim, circular_Ms)
	@test set_Ms(sim, 8.6e5)

	add_exch(sim, 1.3e-11, name="exch")
    add_zeeman(sim, (0,0,4e5))
    add_dmi(sim, 4e-3, name="dmi")

    init_m0(sim, m0_fun)
    relax(sim, maxsteps=1000, stopping_torque = 0.1, save_vtk_every=-1)
end

test_energy_minimization(mesh)
