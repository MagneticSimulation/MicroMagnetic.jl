using JuMag
using Printf
using NPZ

mesh =  FDMeshGPU(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9)

function relax_system()
  sim = Sim(mesh, name="std4_relax")
  set_Ms(sim, 8.0e5)
  sim.driver.alpha = 0.5
  sim.driver.precession = false

  add_exch(sim, 1.3e-11)
  add_demag(sim)

  init_m0(sim, (1, 0.25, 0.1))

  relax(sim, maxsteps=5000, stopping_dmdt = 0.1, save_m_every=1)
  npzwrite("m0.npy", sim.spin)
end

relax_system()
#run_dynamics()
