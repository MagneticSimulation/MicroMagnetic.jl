# ---
# title: Standard Problem 4 (sim_with)
# author: Weiwei Wang
# description: an example for the micromagnetic standard problem 4 using MicroMagnetic.jl using high-level interface.
# link: <https://www.ctcms.nist.gov/~rdm/mumag.org.html/>
# ---

# In this example, we demonstrate how to use **MicroMagnetic.jl** to simulate standard problem 4, which involves a thin film geometry. 
# We will run the simulation on a GPU, process the results, and visualize the magnetization dynamics.

# Import MicroMagnetic.jl
using MicroMagnetic

# Enable GPU acceleration
@using_gpu()



# Define the system geometry: a film with thickness t = 3 nm, length L = 500 nm, and width d = 125 nm.
# Gather all the parameters related to standard problem 4:
args = (
    name = "std4",
    task_s = ["relax", "dynamics"],           # List of tasks
    mesh = FDMesh(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9),
    Ms = 8e5,                                 # Saturation magnetization
    A = 1.3e-11,                              # Exchange constant
    demag = true,                             # Enable demagnetization
    m0 = (1, 0.25, 0.1),                      # Initial magnetization
    alpha = 0.02,                             # Gilbert damping
    steps = 100,                              # Number of steps for dynamics
    dt = 0.01ns,                              # Time step size
    stopping_dmdt = 0.01,                     # Stopping criterion for relaxation
    dynamic_m_interval = 1,                   # Save magnetization at each step
    H_s = [(0,0,0), (-24.6mT, 4.3mT, 0)]      # Static field sweep
);

# Run the simulation using the `sim_with` function:
sim = sim_with(args);

# With the above code, the simulation for standard problem 4 is complete. Next, we proceed to process the data, 
# such as visualizing the magnetization distribution or creating a movie from the simulation results.


# Visualize the magnetization dynamics using CairoMakie
using CairoMakie
# Generate a movie from the simulation results stored in the jld2 file
jld2movie("std4.jld2"; output="assets/std4.mp4", component='x');

# Display the generated movie
# ![](./assets/std4.mp4)


# Finally, we plot the time evolution of the magnetization components to compare the results with OOMMF simulations.


# Import necessary functions to plot the time evolution of magnetization
using DelimitedFiles

function plot_m_ts()
    #Load data
    data = readdlm("std4_llg.txt"; skipstart=2)
    oommf = readdlm("assets/std4_oommf.txt")

    #Create a figure for the plot
    fig = Figure(size=(800, 480))
    ax = Axis(fig[1, 1], xlabel="Time (ns)", ylabel="m")

    #Plot OOMMF results
    lines!(ax, oommf[:, 1] * 1e9, oommf[:, 2], label="OOMMF")
    lines!(ax, oommf[:, 1] * 1e9, oommf[:, 3])
    lines!(ax, oommf[:, 1] * 1e9, oommf[:, 4])

    #Plot MicroMagnetic results
    scatter!(ax, data[:, 2] * 1e9, data[:, 4], markersize=6, label="MicroMagnetic.jl")
    scatter!(ax, data[:, 2] * 1e9, data[:, 5], markersize=6)
    scatter!(ax, data[:, 2] * 1e9, data[:, 6], markersize=6)

    #Add legend to the plot
    axislegend()

    return fig
end

# Plot the magnetization time series
plot_m_ts()