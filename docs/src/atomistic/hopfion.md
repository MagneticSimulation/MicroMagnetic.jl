```@meta
ShareDefaultModule = true
```

# Magnetic hopfion

In this example, we demostrate how to get a magnetic hopfion using MicroMagnetic.jl. We use the parameters given in PRL **124**, 127204 (2020).

````@example
using MicroMagnetic
using CairoMakie

mesh = CubicMesh(nx=64, ny=64, nz=32, dx=1e-9, dy=1e-9, dz=1e-9, pbc="xyz")

function relax_system()
    
    sim = Sim(mesh; driver="SD", name="hopfion")

    set_mu_s(sim, 1)

    hopf = hopfion(R=20e-9, p=1, q=1)
    init_m0(sim, hopf)
    
    J1 = 1
    add_exch(sim, J1; name="exch", J2=-0.164, J3=0, J4=-0.082)

    relax(sim; max_steps=50000, stopping_dmdt=1e-3, using_time_factor=false)
    
    MicroMagnetic.save_vtk_points(sim, "hopfion")

    return sim
end
````

Recall the function `relax_system` to obtain the hopfion.

````@example
sim = relax_system()
nothing #hide
````

After obtain the hopfion, we plot the hopfion

````@example
plot_m(sim, k=16)
````