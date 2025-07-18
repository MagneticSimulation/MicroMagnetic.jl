```@meta
ShareDefaultModule = true
```

# Skyrmion lattice

In this example, we demostrate how to get a skyrmion lattice using MicroMagnetic.jl. We use the parameters given in PRL **108** 017206 (2012), which is dimensionless.

````@example
using MicroMagnetic
using CairoMakie

function m0_fun(i, j, k, dx, dy, dz)
    i0, j0, r = 166, 96, 25
    i1 = i % i0
    j1 = j % j0

    if ((i1 - r)^2 + (j1 - r)^2 < r^2)
        return (0.05, 0.01, -1)
    elseif ((i1 - i0 / 2.0 - r)^2 + (j1 - j0 / 2.0 - r)^2 < r^2)
        return (0.05, 0.01, -1)
    end

    return (0, 0, 1)
end

function relax_system()
    mesh = CubicMesh(; nx=166 * 2, ny=96 * 3, nz=1, pbc="xy")
    sim = Sim(mesh; driver="SD", name="skx_latttice")
    set_mu_s(sim, 1.0)

    #Set the exchange, dmi and zeeman
    add_exch(sim, 1.0; name="exch")
    add_zeeman(sim, (0, 0, 3.75e-3))
    add_dmi(sim, 0.09; name="dmi")

    init_m0(sim, m0_fun)
    relax(sim; max_steps=2000, stopping_dmdt=1e-5)

    save_vtk(sim, "skx_latttice.vts")

    return sim
end
````

Recall the function `relax_system` to obtain the skyrmion lattice.

````@example
sim = relax_system();
````

After obtain the skyrmion, we use the following script to plot the skyrmion

````@example
fig = plot_m(sim)
nothing #hide
````

```@setup
save("../public/skx_lattice.png", fig)
```

![](../public/skx_lattice.png)