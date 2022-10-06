using JuMag

JuMag.cuda_using_double(true)

function fun_dw(i,j,k,dx,dy,dz)
  if i < 150
    return (1,0.1,0)
  elseif i<160
   return (1,1,-0.1)
  else
   return (-1,0.1,0)
  end
end

mesh =  FDMeshGPU(nx=500, ny=20, nz=2, dx=2e-9, dy=2e-9, dz=2e-9)

sim = create_sim(mesh, name="DW", Ms=8e5, A=1.3e-11, Ku=1e5, axis=(1,0,0), m0=fun_dw)
#sim.driver.min_tau = 1e-10

relax(sim, maxsteps=2000, stopping_dmdt=0.03, save_vtk_every = -1)

save_sim_data(sim) 

save_vtk(sim, "dw")
