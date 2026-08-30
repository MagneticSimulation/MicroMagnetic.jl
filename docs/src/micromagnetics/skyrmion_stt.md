```@meta
ShareDefaultModule = true
```

# Skyrmion dynamics under spin transfer torques

In this example, we will study the skyrmion dynamics in a 2d film. We will save the skyrmion
positions to a text file and and generate a movie at the end of simulation.

We import the necessary modules

````@example
using MicroMagnetic
using Printf

@using_gpu()
````

set_verbose()

The studied system is a 800nm x 300nm x 2nm film with periodic boundary conditions.

````@example
mesh = FDMesh(; nx=400, ny=150, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy");
nothing #hide
````

We create a Sim instance with the `Sim` constructor and set basic simulation parameters such as
'Ms', 'A', 'D' and 'H'.

````@example
sim = Sim(mesh; driver="SD", name="skx")
set_Ms(sim, 3.87e5)
add_exch(sim, 5.2e-12)
add_dmi(sim, 1e-3)
add_zeeman(sim, (0, 0, 160 * mT))
nothing #hide
nothing #hide
````

We set the initial magnetization configuration to a skyrmion at position (200nm, 80nm) with radius 20 nm.
Note that the initialized magnetization is a roughly guessing for the skyrimon.

````@example
init_m0(sim, skyrmion(center=(200e-9, 80e-9), R=2e-8))
````

We relax the system to obtain the skyrmion profile.

````@example
relax(sim; max_steps=20000, stopping_dmdt=0.01)
````

We save the magnetization to vtk, which can be opened using Paraview for 3D visualization.

````@example
MicroMagnetic.save_vtk(sim, "skx")
````

After obataining the skyrmion profile, we then move the skyrmion using spin transfer torques.
So we change the driver to "LLG" using the `set_driver` function and add the torque with `add_stt`,
setting the relevant parameters such as alpha and the Zhang-Li model parameters.

````@example
set_driver(sim; driver="LLG", alpha=0.05)
add_stt(sim; model=:zhang_li, b=-20, J=(1, 0, 0), xi=0.2)
````

During the simulation, we need to extract the skyrmion center, so we write a call back function
in which the skyrmion positions are computed using the `compute_guiding_center` function and
saved to a text file with append mode.

````@example
function call_back_fun(sim, t)
    Rx, Ry = compute_guiding_center(sim)
    open("XY.txt", "a") do f
        return write(f, @sprintf("%g  %g  %g\n", t, Rx, Ry))
    end
end
````

We use the `run_sim` function to run the simulation.
The magnetization snapshots are saved as ovf files in the folder `skx_LLG`,
which can be exported to a movie (.mp4) using the `ovf2movie`.

````@example
if !isfile("assets/skx.mp4")
    run_sim(sim; steps=100, dt=1e-10, save_m_every=1, call_back=call_back_fun)
    ovf2movie("skx_LLG"; output="assets/skx.mp4", component='z')
end
````

<video src='../assets/skx.mp4' controls='controls'></video>

