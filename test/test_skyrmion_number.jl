using MicroMagnetic
using Printf
using NPZ
using Test

mesh =  FDMesh(nx=100, ny=100, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")

function m0_fun(i,j,k,dx,dy,dz)
  r2 = (i-50)^2 + (j-50)^2
  if r2 < 10
    return (0.1, 0, -1)
  end
  return (0,0,1)
end

function relax_system()
  sim = Sim(mesh, driver="LLG", name="sim")
  sim.driver.precession = false
  sim.driver.alpha = 0.5
  set_Ms(sim, 8e5)

  add_exch(sim, 1.3e-11, name="exch")
  add_zeeman(sim, (0,0,4e5))
  add_dmi(sim, 4e-3, name="dmi")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=2000, stopping_dmdt=0.1, save_vtk_every = 100, save_m_every=-1)
  #println(sim.spin)
  v = zeros(Float64, sim.n_total)
  compute_skyrmion_number(v, sim.spin, mesh)
  Rxs, Rys = compute_guiding_centre(sim.spin, mesh)
  println(sum(v)," ", Rxs, " ", Rxs)
  return sum(v), Rxs[1], Rys[1]
  #npzwrite("m0.npy", sim.spin)
  #save_vtk(sim, "skx", fields=["exch", "dmi"])
end

Q, Rx, Ry = relax_system()
@test abs(Q+1) < 1e-15
@test abs(Rx-100e-9)<1e-9
@test abs(Ry-100e-9)<1e-9
