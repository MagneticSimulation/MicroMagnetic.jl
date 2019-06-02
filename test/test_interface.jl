using JuMag
using Test

mesh = FDMesh(dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc="xy")

function circular_Ms(i,j,k,dx,dy,dz)
	if i^2 + j^2 <= 10^2
		return 8.6e5
	end
	return 0.0
end

function m0_fun(i,j,k,dx,dy,dz)
  L = 50*dx
  x = i*dx
  return sin(2*pi*x/L), sin(2*pi*x/L+1.2), sin(2*pi*x/L+2.3)
end

function time_fun(t)
    return sin(t)
end

function test_energy_minimization(mesh)
    sim = Sim(mesh, driver="SD", name="sim")
    @test set_Ms(sim, circular_Ms)
    @test set_Ms(sim, 8.6e5)

    add_exch(sim, 1.3e-11, name="exch")
    add_dmi(sim, 4e-3, name="dmi")
    add_anis(sim, 1e2, axis=(0,0,1))
    add_zeeman(sim, (1,2,4e5))
    add_demag(sim)

    init_m0(sim, m0_fun)
    relax(sim, maxsteps=100, stopping_torque = 0.1, save_vtk_every=-1)
end

test_energy_minimization(mesh)
