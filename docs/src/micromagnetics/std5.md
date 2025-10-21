```@meta
ShareDefaultModule = true
```

# Standard Problem 5

link: <https://www.ctcms.nist.gov/~rdm/std5/spec5.xhtml>


This script demonstrates the simulation of a rectangular film of magnetic material
using MicroMagnetic.jl. The film has dimensions 100 nm × 100 nm × 10 nm.
The simulation is divided into two main steps:
1. Relaxing the system to obtain an optimal energy state.
2. Simulating the dynamics of a vortex under an applied current.

````@example
using MicroMagnetic
using CairoMakie
using Printf
using DelimitedFiles

@using_gpu()
````

The system is a rectangular film of magnetic material with dimensions 100 nm × 100 nm × 10 nm.

````@example
mesh = FDMesh(; nx=20, ny=20, nz=2, dx=5e-9, dy=5e-9, dz=5e-9);
nothing #hide
````

Initialize a vortex roughly.

````@example
function init_fun(i, j, k, dx, dy, dz)
    x = i - 10
    y = j - 10
    r = (x^2 + y^2)^0.5
    if r < 2
        return (0, 0, 1)
    end
    return (-y / r, x / r, 0)
end
````

## Step 1: Relax the system.
This function relaxes the system into an energy-optimal state using the steepest descent driver.

````@example
function relax_system(mesh)
    sim = Sim(mesh; driver="SD", name="std5")

    A = 1.3e-11  # Exchange constant
    Ms = 8e5     # Saturation magnetization
    set_Ms(sim, Ms)
    add_exch(sim, A)  # Exchange length=5.7nm
    add_demag(sim)

    init_m0(sim, init_fun)
    relax(sim; max_steps=10000, save_m_every=-1)

    return sim
end

sim = relax_system(mesh);
nothing #hide
````

Plot the magnetization distribution using the plot_m function.

````@example
plot_m(sim; component='z')
````

## Step 2: Vortex Dynamics
We change the driver back to LLG and add spin transfer torques to the system.

````@example
set_driver(sim; driver="LLG", alpha=0.1, gamma=2.211e5)
add_stt(sim, model=:zhang_li, P=1.0, Ms=8e5, xi=0.05, J=(1e12, 0, 0))

# we add a SaverItem to save the guiding center each step
center = SaverItem(("Rx", "Ry"), ("<m>", "<m>"), compute_guiding_center)
run_sim(sim, steps=100, dt=5e-11, save_data=true, saver_item=center)
````

We plot the vortex center as a function of time.

````@example
function plot_center()
    data, unit = read_table("std5_llg.txt")
    fig = Figure(; size=(400, 280), backgroundcolor = :transparent)
    ax = Axis(fig[1, 1]; xlabel="time (ns)", ylabel="Vortex center (nm)", backgroundcolor = :transparent)
    ts = data["time"]*1e9
    Rx, Ry = data["Rx"]*1e9, data["Ry"]*1e9
    scatterlines!(ax, ts, Rx; markersize=6, label="cx")
    scatterlines!(ax, ts, Ry; markersize=6, label="cy")
    axislegend()
    return fig
end

fig = plot_center()
````

```@setup
save("../public/std5_center.png", fig)  #Save the plot
```