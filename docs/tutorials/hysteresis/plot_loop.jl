using DelimitedFiles
using CairoMakie
using MicroMagnetic

function plot_loop()
    data, units = read_table("disk_llg.txt")
    m, H = data["m_x"], data["zeeman_Hx"]/mT
    
    fig = Figure(; size=(400, 280), backgroundcolor = :white)
    ax = Axis(fig[1, 1]; xlabel="H (mT)", ylabel="mx", backgroundcolor = :white)

    scatterlines!(ax, H, m; markersize=8, color=:blue, markercolor=:orange)
    scatterlines!(ax, -H, -m; markersize=8, color=:blue, markercolor=:orange)

    save("loop.png", fig)
    return fig
end

fig = plot_loop()