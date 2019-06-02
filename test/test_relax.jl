using JuMag
using Test
using NPZ

function init_dw(i,j,k,dx,dy,dz)
  if i < 150
    return (1,0.1,0)
  elseif i<160
   return (0,1,0.1)
  else
   return (-1,0.1,0)
  end
end

function relax_system_SD(mesh)
	sim = Sim(mesh, name="relax_sd", driver="SD")
	set_Ms(sim, 8.6e5)

	add_exch(sim, 1.3e-11)
	add_anis(sim, 1e5, axis=(1,0,0))

	init_m0(sim, init_dw)
	relax(sim, maxsteps=2000, stopping_torque=0.1, save_vtk_every = 1000)
end

function relax_system_LLG(mesh)
	sim = Sim(mesh, name="relax_llg", driver="LLG")
    sim.driver.precession = false
    sim.driver.alpha = 0.5
	set_Ms(sim, 8.6e5)

	add_exch(sim, 1.3e-11)
	add_anis(sim, 1e5, axis=(1,0,0))

	init_m0(sim, init_dw)
	relax(sim, maxsteps=2000, stopping_dmdt=0.1, save_vtk_every = 1000)
end


mesh =  FDMesh(nx=500, ny=1, nz=11, dx=2e-9, dy=2e-9, dz=1e-9)
relax_system_SD(mesh)
relax_system_LLG(mesh)
