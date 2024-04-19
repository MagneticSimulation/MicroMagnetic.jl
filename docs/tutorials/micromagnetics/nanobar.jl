# ---
# title: Magnetization state of a nanobar
# author: Weiwei Wang
# date: 2022-10-04
# description: an example to demostrate how to obtain the magnetization distribution in JuMag.
# tag: tutorial
# ---


# In this example, we consider a nanobar with dimensions 60nm x 10nm x 5 nm
# We first import JuMag
using JuMag
using Random
using NPZ

# We create a FDMesh, which means that this program will be run using CPU.
mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);

# We create a Sim wth `SD` driver
sim = Sim(mesh, driver="SD", name="bar");

# We set the saturation magnetization of the system.
set_Ms(sim, 8e5)

# We consider two energies, i.e., exchange and demag 
add_exch(sim, 1.3e-12, name="exch")
add_demag(sim, name="demag");

# We use the random initial state, we set a seed to make our program repeatable
Random.seed!(12345)
init_m0_random(sim);

# Save the initial magnetization state for later comparsion
# npzwrite("bar_init.npy", Array(sim.spin))

# Relax the system to obtain the magnetization distribution. The stopping criteria is the stopping_dmdt, 
# typically its value should be in the rangle of [0.01, 1]. 
relax(sim, maxsteps=2000, stopping_dmdt=0.01) 

# Save the magnetization state for later postprocessing
npzwrite("bar.npy", Array(sim.spin)) 

# Save the vtk as well, which can be visualization using Paraview (https://www.paraview.org/)
save_vtk(sim, "bar", fields=["exch", "demag"])

# We use the following script to plot the magnetization
using CairoMakie

function plot_spatial_m()
    folder = @__DIR__
    m = npzread(folder*"/bar.npy")

    nx, ny, nz = 30, 5, 2
    xs = [i*2 for i=1:nx]
    ys = [j*2 for j=1:ny]
    m = reshape(m, 3, nx, ny, nz)
    mx = m[1,:,:,1]
    my = m[2,:,:,1]
    
    fig = Figure(size = (800, 200))
    ax = Axis(fig[1, 1], backgroundcolor = "white")

    arrows!(ax, xs, ys, mx, my, arrowsize = 10, lengthscale = 2, 
            arrowcolor = vec(my), linecolor = vec(my), align = :center)

    save(folder*"/assets/bar_magnetization.png", fig) #src

    return fig

end

plot_spatial_m()
