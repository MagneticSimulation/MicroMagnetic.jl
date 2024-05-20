# ---
# title: Micromagnetic Standard Problem 4
# author: Weiwei Wang
# date: 2024-05-19
# description: an example for the micromagnetic standard problem 4 using JuMag.
# link: <https://www.ctcms.nist.gov/~rdm/mumag.org.html/>
# ---

# We import the necessary modules
using JuMag
using CairoMakie
using DelimitedFiles

# make use of available gpu 
@using_gpu()

# The geometry of the studied system is a film of thickness t=3 nm, length L=500 nm and width d=125 nm.
mesh =  FDMesh(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);

# The first step is to relax the system, after which we should obtain a 'S' state.
sim = create_sim(mesh, name="std4", driver="SD", Ms=8.0e5, A=1.3e-11, demag=true,  m0=(1, 0.25, 0.1));
relax(sim, maxsteps=5000, stopping_dmdt=0.01)

# We plot the magnetization distribution using the plot_m function
plot_m(sim)

# The second step is to apply an external field starting from the obtained 'S' state in previous step. 
set_driver(sim, driver="LLG", alpha=0.02, gamma = 2.211e5)
add_zeeman(sim, (-24.6mT, 4.3mT, 0))
run_sim(sim, steps=100, dt=1e-11)

# By default, run_sim will generate a jdl2 file as well, so we can generate a movie based on it.
jdl2movie("std4.jdl2", output="assets/std4.mp4")

# ![](./assets/std4.mp4)

# Finally, we plot the time evolution of the magnetization using the plot_m_ts function
function plot_m_ts()
    folder = @__DIR__
    data = readdlm("std4_llg.txt", skipstart=2)
    oommf = readdlm("assets/std4_oommf.txt")

    fig = Figure(size = (800, 480))
    ax = Axis(fig[1, 1],
        xlabel = "Time (ns)",
        ylabel = "m"
    )

    lines!(ax, oommf[:,1]*1e9, oommf[:,2], label="OOMMF")
    lines!(ax, oommf[:,1]*1e9, oommf[:,3])
    lines!(ax, oommf[:,1]*1e9, oommf[:,4])

    scatter!(ax, data[:,2]*1e9, data[:,4], markersize = 6, label="JuMag")
    scatter!(ax, data[:,2]*1e9, data[:,5], markersize = 6)
    scatter!(ax, data[:,2]*1e9, data[:,6], markersize = 6)

    axislegend()

    return fig

end

plot_m_ts()