```@meta
ShareDefaultModule = true
```

# Standard Problem 4

This example demonstrates how to simulate standard problem 4 using **MicroMagnetic.jl**. We start by relaxing the system to achieve a stable magnetization configuration,
followed by applying an external magnetic field to study its effects.

Import necessary modules

````@example
using MicroMagnetic
using CairoMakie
using DelimitedFiles
````

Enable GPU acceleration

````@example
@using_gpu()
````

Define the system geometry: a film with thickness t = 3 nm, length L = 500 nm, and width d = 125 nm.

````@example
mesh = FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);
nothing #hide
````

## Step 1: Relaxing the System
The first step in our simulation is to relax the system to obtain a stable magnetization configuration, often referred to as an "S" state. We encapsulate this process in the `relax_system` function.

````@example
function relax_system(mesh)
   sim = Sim(mesh; driver="SD", name="std4")

   Ms = 8e5
   A = 1.3e-11

   set_Ms(sim, Ms)   # Set saturation magnetization
   add_exch(sim, A)  # Add exchange interaction
   add_demag(sim)    # Add demagnetization

   init_m0(sim, (1, 0.25, 0.1))  # Initialize magnetization
   relax(sim; stopping_dmdt=0.01)  # Relax the system

   return sim
end
````

The `relax_system` function takes a mesh as input and performs the following steps:

**Simulation Initialization:**
The function initializes a simulation (`sim`) using the given mesh and sets the driver to "SD" (Steepest Descent), which is typically used for relaxation processes. The simulation is named "std4" for consistency with standard problem 4.

**Material Parameters Setup:**
The saturation magnetization `Ms` and exchange constant `A` are set up using `set_Ms` and `add_exch` functions, respectively. Demagnetization effects are added using `add_demag`.

**Initial Magnetization:**
The initial magnetization vector `m0` is set using `init_m0`. In this case, it is initialized with a vector `(1, 0.25, 0.1)`.

**Relaxation:**
The system is relaxed using the `relax` function, which iteratively minimizes the system's energy until the change in magnetization (`dm/dt`) falls below a specified threshold (`stopping_dmdt=0.01`).

Relax the system to obtain a stable magnetization configuration

````@example
sim = relax_system(mesh);
nothing #hide
````

## Step 2: Applying an External Field
After obtaining the stable "S" state, the next step is to apply an external magnetic field and observe the magnetization dynamics.
Plot the magnetization distribution using the `plot_m` function

````@example
plot_m(sim; component='x')
````

Apply an external field starting from the relaxed "S" state

```julia
set_driver(sim; driver="LLG", alpha=0.02, gamma=2.211e5)
add_zeeman(sim, (-24.6mT, 4.3mT, 0))  # Apply external magnetic field
run_sim(sim; steps=100, dt=1e-11, save_m_every=1)   # Run the simulation for 10 steps
```

## Step 3: Visualizing the Results
By default, `run_sim` generates a `std4_LLG` folder, which stores the magnetization data. You can create a movie from it to visualize the magnetization dynamics.
Generate a movie based on the simulation results

```julia
ovf2movie("std4_LLG"; output="../public/std4.gif", component='x');
```
![](../public/std4.gif)