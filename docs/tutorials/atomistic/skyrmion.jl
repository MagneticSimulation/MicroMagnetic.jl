# ---
# title: Magnetic skyrmion
# author: Weiwei Wang
# date: 2022-10-04
# description: an example to demostrate how to obtain a skyrmion in JuMag.
# tag: atomistic; skyrmion
# ---

# In this example, we show how to obtain a magnetic skyrmion using atomistic module in JuMag.
# We import JuMag and use double float precision in the simulation.
using JuMag
using CairoMakie

# We create a CubicMesh
mesh =  CubicMesh(nx=120, ny=120, nz=1, pbc="xy");

# To start the simulation, we need to give an initial state.
# We define a function in which we set the spins around site (60,60) to be negative
function m0_fun(i,j,k, dx, dy, dz)
  r2 = (i-60)^2 + (j-60)^2
  if r2 < 10^2
    return (0.1, 0, -1)
  end
  return (0,0,1)
end


# We define a function to specify the problem.
function relax_system(mesh)
  #We create a simulation with 'SD' driver
  sim = Sim(mesh, driver="SD", name="skyrmion")

  #We set mu_s of the system
  set_mu_s(sim, 1.0)
  sim.driver.max_tau = 1000.0

  #Set the exchange, dmi and zeeman
  add_exch(sim, 1.0, name="exch")
  add_zeeman(sim, (0,0,3.75e-3))
  add_dmi(sim, 0.09, name="dmi")

  #Initialize the system using the `m0_fun` function
  init_m0(sim, m0_fun)

  #Relax the system
  relax(sim, maxsteps=2000, stopping_dmdt=1e-4, using_time_factor=false)

  return sim
end

# Recall the function 
sim = relax_system(mesh);

# After obtain the skyrmion, we use the following script to plot the skyrmion
plot_m(sim)


