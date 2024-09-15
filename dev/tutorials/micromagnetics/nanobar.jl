using MicroMagnetic
using CairoMakie

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);

sim = Sim(mesh; driver="SD")
set_Ms(sim, 8e5)   #Set saturation magnetization Ms=8e5 A/m

add_exch(sim, 1e-12);  #Add exchange interaction
add_demag(sim);         #Add demagnetization

init_m0(sim, (1, 1, 0))  #Initialize magnetization

plot_m(sim)

relax(sim; max_steps=2000, stopping_dmdt=0.01)

fig = plot_m(sim)

save("bar.png", fig)

save_vtk(sim, "bar"; fields=["exch", "demag"])

args = (
    task = "Relax",
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2),
    Ms = 8e5,
    A = 1e-12,
    demag = true,
    m0 = (1, 1, 0),
    stopping_dmdt = 0.01
);

sim = sim_with(args);

plot_m(sim)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
