using JuMag

JuMag.cuda_using_double(true)

function init_dw(i,j,k,dx,dy,dz)
  if i < 150
    return (1,0.1,0)
  elseif i<160
   return (0,1,0.1)
  else
   return (-1,0.1,0)
  end
end

mesh =  FDMeshGPU(nx=500, ny=2, nz=2, dx=2e-9, dy=4e-9, dz=4e-9)

sim = StdSim(mesh, Ms=8e5, A=1.3e-11, Ku=1e5, axis=(1,0,0), name="DW", m0=init_dw)

relax(sim, maxsteps=2000, stopping_dmdt=0.03, save_vtk_every = -1)

save_vtk(sim, "dw")
