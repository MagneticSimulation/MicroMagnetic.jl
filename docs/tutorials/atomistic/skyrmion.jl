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
using NPZ
JuMag.cuda_using_double(true)

# We create a CubicMesh
mesh =  CubicMeshGPU(nx=120, ny=120, nz=1, pbc="xy")

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

  #Set the max step 
  sim.driver.max_tau = 1.0

  #Set the exchange, dmi and zeeman
  add_exch(sim, 1.0, name="exch")
  add_zeeman(sim, (0,0,3.75e-3))
  add_dmi(sim, 0.09, name="dmi")

  #Initialize the system using the `m0_fun` function
  init_m0(sim, m0_fun)

  #Relax the system
  relax(sim, maxsteps=2000, stopping_dmdt=200)

  #Save the final magnetization state for later postprocessing
  npzwrite("skx.npy", Array(sim.spin))
  
  #Save the vtk as well
  save_vtk(sim, "skx", fields=["exch", "dmi"])
end

# Recall the function 
relax_system(mesh)

# After obtain the skyrmion, we use the following script to plot the skyrmion
using CairoMakie

function plot_spatial_m()
  m = npzread("skx.npy")

  nx, ny, nz = 120, 120, 1
  points = [Point3f(i, j, 0) for i in 1:5:nx for j in 1:5:ny]
  
  m = reshape(m, 3, nx, ny)
  mf = [Vec3f(m[1, i, j], m[2, i,j], m[3, i,j]) for i in 1:5:nx for j in 1:5:ny]
  mz = [m[3, i, j]  for i in 1:5:nx for j in 1:5:ny]

  fig = Figure(resolution = (1000, 1000))
  ax = Axis(fig[1, 1], backgroundcolor = "white")

  arrows!(ax, points, mf, fxaa=true, # turn on anti-aliasing
          color = vec(mz), linewidth = 1, arrowsize = 2, lengthscale = 2,
          align = :center
      )

  save("assets/skx.png", fig, px_per_unit = 1) #src

  return fig

end

plot_spatial_m()


