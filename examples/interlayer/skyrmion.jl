using JuMag
using Printf
using NPZ

mesh =  FDMeshGPU(nx=100, ny=100, nz=3, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")

function m0_fun(i,j,k,dx,dy,dz)
  r2 = (i-50)^2 + (j-50)^2
  if r2 < 100
    return (0.1, 0, -1)
  end
  return (0,0,1)
end

function spatial_Ms(i,j,k,dx,dy,dz)
  if k == 2
    return 0
  else
    return 8e5
  end
end

function relax_system()
  sim = Sim(mesh, driver="SD", name="sim")
  set_Ms(sim, spatial_Ms)

  add_exch(sim, 1.3e-11, name="exch")
  add_zeeman(sim, (0,0,4e5))
  add_dmi(sim, 4e-3, name="dmi")
  add_dmi_interlayer(sim, (0.0, 0.0, 1e-4), name="interlayer")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=1000, stopping_dmdt=0.1, save_vtk_every = 10)

  save_vtk(sim, "skx", fields=["exch", "dmi", "interlayer"])
end

relax_system()

