using JuMag
using Printf
using NPZ

mesh =  CubicMesh(nx=100, ny=100, nz=1, pbc="xy")

function m0_fun(i,j,k)
  r2 = (i-50)^2 + (j-50)^2
  if r2 < 10^2
    return (0.1, 0, -1)
  end
  return (0,0,1)
end

function relax_system()
  sim = Sim(mesh, driver="LLG", name="sim")
  set_mu_s(sim, 1.0)
  sim.driver.gamma = 1.0
  sim.driver.alpha = 0.5
  sim.driver.precession = false

  add_exch(sim, 1.0, name="exch")
  add_zeeman(sim, (0,0,3.75e-3))
  add_dmi(sim, 0.09, name="dmi")

  init_m0(sim, m0_fun)

  relax(sim, maxsteps=2000, stopping_dmdt=1e-4, save_vtk_every = 100)
  #println(sim.spin)
  npzwrite("m0.npy", sim.spin)
  #save_vtk(sim, "skx", fields=["exch", "dmi"])
end

function relax_system_sdm()
  sim = Sim(mesh, driver="SD", name="sim_sd")
  set_mu_s(sim, 1.0)
  sim.driver.max_tau = 1.0

  add_exch(sim, 1.0, name="exch")
  add_zeeman(sim, (0,0,3.75e-3))
  add_dmi(sim, 0.09, name="dmi")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=2000, stopping_dmdt=0.01, save_vtk_every = 10)
  #println(sim.spin)
  npzwrite("m2.npy", sim.spin)
  #save_vtk(sim, "skx", fields=["exch", "dmi"])
end

relax_system()
