```@meta
ShareDefaultModule = true
```

# Skyrmion lattice

In this example, we demostrate how to get a magnetic hopfion using MicroMagnetic.jl. We use the parameters given in PRL **124**, 127204 (2020).

````@example
using MicroMagnetic
using CairoMakie

mesh = CubicMesh(nx=64, ny=64, nz=32, dx=1e-9, dy=1e-9, dz=1e-9)

function init_hopfion(x, y, z)
    r = sqrt(x^2+y^2)
    d = 10e-9
    sinf = 2*r*d/(r^2+d^2)
    cosf = sqrt(1-sinf^2)
    mx = x/r*2*sinf*cosf + y*z/r^2*sinf^2
    my = y/r*2*sinf*cosf - x*z/r^2*sinf^2
    mz = cosf^2-sinf^2 + 2*z^2/r^2*sinf^2
    return (mx, my, mz)
end

function relax_system()
    
    sim = Sim(mesh; driver="SD", name="hopfion")

    set_mu_s(sim, 1)

    init_m0(sim, init_hopfion)

    J1 = 1
    add_exch(sim, J1; name="exch", J2=-0.164, J3=0, J4=-0.082)

    relax(sim; max_steps=50000, stopping_dmdt=5e-3, using_time_factor=false)
    
    MicroMagnetic.save_vtk_points(sim, "hopfion")
end
````

Recall the function `relax_system` to obtain the hopfion.

````@example
sim = relax_system()
nothing #hide
````

After obtain the hopfion, we use the following script to plot the hopfion

````@example
fig = plot_m(sim, k=32)
````