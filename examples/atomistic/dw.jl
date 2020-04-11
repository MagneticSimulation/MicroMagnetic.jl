using JuMag
using Printf
using Test
using NPZ

function init_dw(i,j,k, dx, dy, dz)
  if i < 150
    return (1,0.1,0)
  elseif i<160
   return (0,1,0.1)
  else
   return (-1,0.1,0)
  end
end

mesh =  CubicMeshGPU(nx=300)
sim = Sim(mesh, name="dw_atomic")

set_mu_s(sim, 1.0)
sim.driver.alpha = 0.5
sim.driver.gamma = 1.0
sim.driver.precession = false

add_exch(sim, 1.0)
add_anis(sim, 0.005, axis=(1,0,0))

init_m0(sim, init_dw)

relax(sim, stopping_dmdt=1e-5, maxsteps=1000)

npzwrite("dw.npy", Array(sim.spin))
