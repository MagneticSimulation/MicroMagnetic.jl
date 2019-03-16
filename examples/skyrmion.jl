using SpinDynamics
using Printf
using NPZ

mesh =  create_mesh(nx=100, ny=100, nz=100, dx=2.0, dy=2.0, dz=2.0, unit_length=1e-9, pbc="xy")

function m0_fun(i,j,k,dx,dy,dz)
  r2 = (i-50)^2 + (j-50)^2
  if r2 < 100
    return (0.0001, 0, -1)
  end
  return (0,0,1)
end

function relax_system()
  sim = create_sim(mesh, name="skx", tol=1e-6)
  sim.Ms[:] .= 8.6e5
  sim.alpha = 0.5
  sim.gamma = 2.211e5
  add_exch(sim, 1.3e-11, name="exch")
  add_zeeman(sim, (0,0,4e5))
  add_dmi(sim, 4e-3, name="dmi")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=5000, stopping_dmdt = 0.1)
  #println(sim.spin)
  npzwrite("m0.npy", sim.spin)
  save_vtk(sim, "skx", fields=["exch", "dmi"])
end

relax_system()
