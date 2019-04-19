using JuMag
using Printf
using NPZ

mesh =  FDMeshGPU(nx=100, ny=100, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")

function m0_fun(i,j,k,dx,dy,dz)
  r2 = (i-50)^2 + (j-50)^2
  if r2 < 100
    return (0.1, 0, -1)
  end
  return (0,0,1)
end

function relax_system()
  sim = Sim(mesh, driver="SD", name="sim")
  set_Ms(sim, 8e5)

  add_exch(sim, 1.3e-11, name="exch")
  add_zeeman(sim, (0,0,4e5))
  add_dmi(sim, 4e-3, name="dmi")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=200, stopping_torque=0.1, save_vtk_every = 10)
  #println(sim.spin)
  npzwrite("m0.npy", sim.spin)
  #save_vtk(sim, "skx", fields=["exch", "dmi"])
end

relax_system()
