using SpinDynamics
using Printf
using Test

function init_dw(i,j,k,dx,dy,dz)
  if i < 150
    return (1,0.1,0)
  elseif i<160
   return (0,1,0.1)
  else
   return (-1,0.1,0)
  end
end

mesh =  create_mesh(nx=300, dx=2e-9)
sim = create_sim(mesh, name="dyn")

sim.Ms[:] .= 8.6e5
sim.alpha = 0.5
sim.gamma = 2.21e5

add_exch(sim, 1.3e-11)
add_anis(sim, 1e5, axis=(1,0,0))

init_m0(sim, init_dw)

relax(sim, maxsteps=500)

using Plots
gr()

xs = 1:300
b = reshape(sim.spin, 3, sim.nxyz)
plot(xs, b[1,:], marker=:h, markersize=3, linewidth=1, label="mx")
plot!(xs, b[2,:], marker=:c, markersize=3, linewidth=1, label="my")
plot!(xs, b[3,:], marker=:p, markersize=1.2, linewidth=1, label="mz")
savefig("dw.png")
