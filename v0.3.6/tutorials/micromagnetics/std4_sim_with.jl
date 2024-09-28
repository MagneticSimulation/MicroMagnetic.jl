using MicroMagnetic

@using_gpu()

args = (name="std4", task_s=["relax", "dynamics"],           # List of tasks
        mesh=FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9), Ms=8e5,                                 # Saturation magnetization
        A=1.3e-11,                              # Exchange constant
        demag=true,                             # Enable demagnetization
        m0=(1, 0.25, 0.1),                      # Initial magnetization
        alpha=0.02,                             # Gilbert damping
        steps=100,                              # Number of steps for dynamics
        dt=0.01ns,                              # Time step size
        stopping_dmdt=0.01,                     # Stopping criterion for relaxation
        dynamic_m_interval=1,                   # Save magnetization at each step
        H_s=[(0, 0, 0), (-24.6mT, 4.3mT, 0)]);

sim = sim_with(args);

using CairoMakie

jld2movie("std4.jld2"; output="assets/std4.mp4", component='x');

using DelimitedFiles

function plot_m_ts()
    #Load data
    data, unit = read_table("std4_llg.txt")
    oommf = readdlm("assets/std4_oommf.txt")

    #Create a figure for the plot
    fig = Figure(; size=(800, 480))
    ax = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="m")

    #Plot OOMMF results
    lines!(ax, oommf[:, 1] * 1e9, oommf[:, 2]; label="OOMMF")
    lines!(ax, oommf[:, 1] * 1e9, oommf[:, 3])
    lines!(ax, oommf[:, 1] * 1e9, oommf[:, 4])

    #Plot MicroMagnetic results
    scatter!(ax, data["time"] * 1e9, data["m_x"]; markersize=6, label="MicroMagnetic.jl")
    scatter!(ax, data["time"] * 1e9, data["m_y"]; markersize=6)
    scatter!(ax, data["time"] * 1e9, data["m_z"]; markersize=6)

    #Add legend to the plot
    axislegend()

    return fig
end

plot_m_ts()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
