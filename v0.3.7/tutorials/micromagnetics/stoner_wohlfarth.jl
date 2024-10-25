using MicroMagnetic

mesh = FDMesh(; nx=4, ny=4, nz=4, dx=1e-9, dy=1e-9, dz=1e-9);

args = (
    name = "sw",
    task = "Relax",
    mesh = mesh,
    Ms=1.0e6,
    A=1.3e-11,
    m0=(-1, 1, 0),
    Ku=5e4,
    axis=(1, 1, 0),
    stopping_dmdt = 0.05,
    H_s = [(i*1mT, 0, 0) for i=-100:5:100]
);

sim_with(args);

using DelimitedFiles
using CairoMakie

function plot_loop()
    data = readdlm("./sw_sd.txt"; skipstart=2)
    m, H = data[:, 3], data[:, 8]

    fig = Figure(; size=(600, 400))
    ax = Axis(fig[1, 1]; xlabel="H (A/m)", ylabel="mx")

    scatterlines!(ax, H, m; markersize=8, color=:blue, markercolor=:orange)
    scatterlines!(ax, -H, -m; markersize=8, color=:blue, markercolor=:orange)

    expected = 39788.736 # A/m
    vlines!(ax, [expected, -expected]; color=:red, linestyle=:dash)

    return fig
end

plot_loop()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
