# ---
# title: Magnetization state of a nanobar
# author: Weiwei Wang
# date: 2024-05-16
# description: an example to demostrate how to obtain the magnetization distribution in JuMag.
# tag: tutorial
# ---


# In this example, we consider a nanobar with dimensions 60nm x 10nm x 5 nm
# We first import JuMag and CairoMakie for plotting.
using JuMag
using CairoMakie

# We create a FDMesh
mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);

# We create a Sim wth `SD` driver using `create_sim` function.
# We consider two energies (i.e., exchange and demag) with exchange constant A = 1e-12 J/m.
# The saturation magnetization Ms=8e5 A/m
sim = create_sim(mesh, driver="SD", name="bar", Ms=8e5, A=1e-12, demag=true, m0=(1, 1, 0));

# We can plot the magnetization using `plot_m` function
plot_m(sim)

# We relax the system to obtain the magnetization distribution. The stopping criteria is the stopping_dmdt, 
# typically its value should be in the rangle of [0.01, 1]. 
relax(sim, maxsteps=2000, stopping_dmdt=0.01) 

# We plot the magnetization
fig = plot_m(sim)

# We save the figure to png. 
save("bar.png", fig)

# Save the magnetization state for later postprocessing, which can be visualization using Paraview (https://www.paraview.org/)
save_vtk(sim, "bar", fields=["exch", "demag"])
