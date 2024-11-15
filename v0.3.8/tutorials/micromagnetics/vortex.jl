using MicroMagnetic
using CairoMakie

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=5e-9, nx=100, ny=100, nz=4);

geo = Cylinder(; radius=100e-9);

function init_fun(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end

args = (
    task = "Relax",
    mesh = mesh,
    shape = geo,
    Ms = 8e5,
    A = 1.3e-11,
    demag = true,
    m0 = init_fun,
    stopping_dmdt = 0.01
);

sim = sim_with(args);

plot_m(sim; figsize=(400, 400), arrows=(30, 30))

save_vtk(sim, "vortex")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
