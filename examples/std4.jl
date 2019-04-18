using JuMag
using Printf
using NPZ

mesh =  FDMeshGPU(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9)

function relax_system(mesh)
  sim = Sim(mesh, name="std4_relax")
  set_Ms(sim, 8.0e5)
  sim.driver.alpha = 0.5
  sim.driver.precession = false

  add_exch(sim, 1.3e-11)
  add_demag(sim)

  init_m0(sim, (1, 0.25, 0.1))

  relax(sim, maxsteps=5000, stopping_dmdt = 0.1, save_m_every=1)
  npzwrite("m0.npy", Array(sim.spin))
end

function apply_field1(mesh)
  sim = Sim(mesh, name="std4_dyn")
  set_Ms(sim, 8.0e5)
  sim.driver.alpha = 0.02
  sim.driver.gamma = 2.211e5

  mT = 0.001 / (4*pi*1e-7)
  add_exch(sim, 1.3e-11)
  add_demag(sim)
  add_zeeman(sim, (-24.6*mT, 4.3*mT, 0))

  init_m0(sim, npzread("m0.npy"))

  for i=1:100
    run_until(sim, 1e-11*i)
  end
end

#relax_system(mesh)
println("Start step2 !!!")
apply_field1(mesh)
#println("Run step2 again!!!")
#@time apply_field1(mesh)
println("Done!")
