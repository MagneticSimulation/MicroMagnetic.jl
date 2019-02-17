using SpinDynamics
using Printf
using NPZ

mesh =  create_mesh(nx=200, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)

function relax_system()
  sim = create_sim(mesh, name="relax", tol=1e-6)
	sim.Ms[:] .= 8.0e5
	sim.alpha = 1.0
	sim.gamma = 2.211e5
  add_exch(sim, 1.3e-11)
  add_demag(sim)

  init_m0(sim, (1, 0.25, 0.1))

  relax(sim, maxsteps=5000, stopping_dmdt = 0.01, save_m_every=1)
  npzwrite("m0.npy", sim.spin)
end

function run_dynamics()
  sim = create_sim(mesh, name="std4", tol=1e-7)
	sim.Ms[:] .= 8.0e5
	sim.alpha = 0.02
	sim.gamma = 2.211e5

  m0 = npzread("m0.npy")
  init_m0(sim, m0)
	add_exch(sim, 1.3e-11)
  add_demag(sim)

  mT = 0.001 / (4*pi*1e-7)
  add_zeeman(sim, (-24.6 * mT, 4.3 * mT, 0))

  for i = 1:100
    run_until(sim, 1e-11*i)
    println("step=",i)
  end
end


relax_system()
#run_dynamics()
