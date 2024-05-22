# ---
# title: Skyrmion lattice
# author: Weiwei Wang
# date: 2022-10-07
# description: to obtain the skyrmion lattice using the parameters shown in PRL **108** 017206 (2012).
# tag: atomistic; skyrmion
# ---

using MicroMagnetic
using NPZ

function m0_fun(i, j, k, dx, dy, dz)
  i0, j0, r = 166, 96, 25
  i1 = i % i0
  j1 = j % j0

  if ((i1 - r)^2 + (j1 - r)^2 < r^2)
    return (0.05, 0.01, -1)
  elseif ((i1 - i0 / 2. - r)^2 + (j1 - j0 / 2. - r)^2 < r^2)
    return (0.05, 0.01, -1)
  end
  
  return (0,0,1)
end

function relax_system()
  mesh =  CubicMesh(nx=166*2, ny=96*3, nz=1, pbc="xy")
  sim = Sim(mesh, driver="SD", name="skx")
  set_mu_s(sim, 1.0)

  add_exch(sim, 1.0, name="exch")
  add_zeeman(sim, (0,0,3.75e-3))
  add_dmi(sim, 0.09, name="dmi")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=2000, stopping_dmdt=0.2, save_vtk_every = -1)

  npzwrite("skx.npy", Array(sim.spin))
  save_vtk(sim, "skx")
end

relax_system()

# Plot the obtained spin structure.
using CairoMakie

function plot_spatial_m()
  m = npzread("skx.npy")

  nx, ny, nz = 166*2, 96*3, 1
  points = [Point3f(i, j, 0) for i in 1:5:nx for j in 1:5:ny]
  
  m = reshape(m, 3, nx, ny)
  mf = [Vec3f(m[1, i, j], m[2, i,j], m[3, i,j]) for i in 1:5:nx for j in 1:5:ny]
  mz = [m[3, i, j]  for i in 1:5:nx for j in 1:5:ny]

  fig = Figure(resolution = (1600, 1600))
  ax = Axis(fig[1, 1], backgroundcolor = "white")

  arrows!(ax, points, mf, fxaa=true, # turn on anti-aliasing
          color = vec(mz), linewidth = 1, arrowsize = 2, lengthscale = 2,
          align = :center
      )

  #save("assets/skx.png", fig, px_per_unit = 1) #src

  return fig

end

plot_spatial_m()
