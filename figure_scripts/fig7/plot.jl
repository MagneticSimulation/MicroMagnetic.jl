using CairoMakie
using DelimitedFiles

function plot_m_H()

  fig = Figure(size = (500, 300))
  
  ax1 = Axis(fig[1, 1],
      xlabel = "T (K)",
      ylabel = "m"
  )

  ax2 = Axis(fig[1, 1],
      xlabel = "T (K)",
      ylabel = L"$\chi$",
      yaxisposition = :right
  )

  ax2.yticklabelcolor = :blue
  ax2.ylabelcolor = :blue

  data = readdlm("assets/M_T.txt", skipstart=1)
  
  sc1 = scatter!(ax1, data[:, 1], data[:, 2], marker=:rect, markersize = 10)
  sc1.color = :transparent
  sc1.strokewidth = 1
  sc1.strokecolor = :purple
  lines!(ax1, data[:, 1], data[:, 2], color=:purple)
  

  sc2 = scatter!(ax2, data[:, 1], data[:, 3], markersize = 8)
  sc2.color = :transparent
  sc2.strokewidth = 1
  sc2.strokecolor = :blue
  lines!(ax2, data[:, 1], data[:, 3], color=:blue)


  vlines!(ax1, 431, linewidth=1, color=:black, linestyle=:dash)
  
  axislegend(ax1, [sc1, sc2], ["m", L"$\chi$"], position=:lc)
  
  save("M_T.pdf", fig)
  save("M_T.png", fig)

  return fig

end

plot_m_H()
