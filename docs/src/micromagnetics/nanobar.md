```@meta
ShareDefaultModule = true
```

# Nanobar

In this example, we consider a nanobar with dimensions 60nm x 10nm x 5 nm
We used this example to demostrate how to obtain the magnetization distribution using MicroMagnetic.jl
We first import MicroMagnetic and CairoMakie for plotting.

````@example
using MicroMagnetic
using CairoMakie
````

We create a FDMesh

````@example
mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);
nothing #hide
````

We create a Sim wth `SD` driver using `Sim` function and set the saturation magnetization Ms

````@example
sim = Sim(mesh; driver="SD")
set_Ms(sim, 8e5)   #Set saturation magnetization Ms=8e5 A/m
````

We consider two energies (i.e., exchange and demag) with exchange constant A = 1e-12 J/m.

````@example
add_exch(sim, 1e-12);   #Add exchange interaction
add_demag(sim);         #Add demagnetization
nothing #hide
````

Initilize the magnetization to (1,1,0) direction,

````@example
init_m0(sim, (1, 1, 0))  #Initialize magnetization
````

We can plot the magnetization using `plot_m` function

````@example
plot_m(sim)
````

We relax the system to obtain the magnetization distribution. The stopping criteria is the stopping_dmdt,
typically its value should be in the rangle of [0.01, 1].

````@example
relax(sim; max_steps=2000, stopping_dmdt=0.01)
````

We plot the magnetization again

````@example
fig = plot_m(sim)
````

We save the figure to png.

```julia
save("bar.png", fig)
```

Save the magnetization state for later postprocessing, which can be visualization using Paraview (https://www.paraview.org/)

```julia
save_vtk(sim, "bar"; fields=["exch", "demag"])
```

## Using the sim_with function.
We can use sim_with to simplify the setup of the simulation.
We put all the parameters together:

````@example
args = (
    task = "Relax",
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2),
    Ms = 8e5,
    A = 1e-12,
    demag = true,
    m0 = (1, 1, 0),
    stopping_dmdt = 0.01
);
nothing #hide
````

then use the sim_with function

````@example
sim = sim_with(args);
nothing #hide
````

We plot the magnetization

````@example
plot_m(sim)
````
