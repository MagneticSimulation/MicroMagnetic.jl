using CairoMakie
# Import necessary functions to plot the time evolution of magnetization
using DelimitedFiles

function plot_m_ts()
    #Load data
    data = readdlm("std4_llg.txt"; skipstart=2)
    oommf = readdlm("../../docs/tutorials/micromagnetics/assets/std4_oommf.txt")

    #Create a figure for the plot
    fig = Figure(size = (500, 360), fontsize = 18)
    ax = Axis(fig[1, 1], xlabel="Time (ns)", ylabel="m")

    #Plot OOMMF results
    l1 = lines!(ax, oommf[:, 1] * 1e9, oommf[:, 2])
    l2 = lines!(ax, oommf[:, 1] * 1e9, oommf[:, 3], color=:sienna1)
    l3 = lines!(ax, oommf[:, 1] * 1e9, oommf[:, 4], color=:slateblue1)

    #Plot MicroMagnetic results
    scatter!(ax, data[:, 2] * 1e9, data[:, 4], markersize=6, label="MicroMagnetic.jl")
    scatter!(ax, data[:, 2] * 1e9, data[:, 5], markersize=6, color=:sienna1)
    scatter!(ax, data[:, 2] * 1e9, data[:, 6], markersize=5, color=:slateblue1)

    #Add legend to the plot
    axislegend(ax, [l1, l2, l3], [L"m_x", L"m_y", L"m_z"], position = :rt, orientation = :horizontal, labelsize=18)

    save("std4.pdf", fig)
    save("std4.png", fig)

    return fig
end

# Plot the magnetization time series
plot_m_ts()