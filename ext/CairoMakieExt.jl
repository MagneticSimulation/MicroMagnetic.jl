
module CairoMakieExt

using JuMag
using Printf
using CairoMakie
using JLD2

#TODO: find a better sampling method
function calculate_sampling(nx::Int, step::Int)
    start = step > 1 ? (step ÷ 2) : 1
    num_samples = (nx - start) ÷ step + 1
    return start, step, num_samples
end

function calculate_start_step(nx::Int, n::Int)
    step = nx > n ? div(nx, n) : 1
    start = (nx - step * (n - 1)) ÷ 2 + 1
    if start <= 0
        start = 1
    end
    return start, step
end

"""
    plot_m(spin; dx=1.0, dy=1.0, k=1, arrows=(-1, 1), figsize=500, fig=nothing)

Create a plotting for given magnetization. 

  `spin` should be an array with dimension (3, nx, ny, nz)
  `k` indicates the layer index (starting from 1) 
  `arrows` is the number of arrows, should be a Tuple of integers. By default, arrows=(-1, -1).
  `figsize` should be a Tuple of integers, for example, figsize=(500, 400) or figsize=(500, -1).
"""
function JuMag.plot_m(spin; dx=1.0, dy=1.0, k=1, arrows=(-1, -1), figsize=(500, -1), fig=nothing, ax=nothing, colorrange=nothing)
    (_, nx, ny, nz) = size(spin)
    scale_factor = 10^floor(log10(dx))
    dx = dx / scale_factor
    dy = dy / scale_factor
    xs = [i * dx for i in 1:nx]
    ys = [j * dy for j in 1:ny]

    mx = spin[1, :, :, k]
    my = spin[2, :, :, k]
    mz = spin[3, :, :, k]
    lml = sqrt.(mx.^ 2 .+ my.^ 2 .+ mz.^2)
    mx[lml .< 0.1] .= NaN
    my[lml .< 0.1] .= NaN
    mz[lml .< 0.1] .= NaN

    size_x = figsize[1]
    size_y = figsize[2]
    if (size_y < 0)
        aspect_ratio = ny * dy / (nx * dx)
        size_y = Int(ceil(size_x * aspect_ratio))
    end

    max_arrows = 40

    arrow_nx = arrows[1]
    arrow_ny = arrows[2]

    if arrows[1] < 0 && arrows[2] > 0
        start_y, step_y = calculate_start_step(ny, arrow_ny)
        start_x, step_x, arrow_nx = calculate_sampling(nx, step_y)
    elseif arrows[1] > 0 && arrows[2] < 0
        start_x, step_x = calculate_start_step(nx, arrow_nx)
        start_y, step_y, arrow_ny = calculate_sampling(ny, step_x)
    elseif arrows[1] < 0 && arrows[2] < 0
        _, step_size = calculate_start_step(max(nx, ny), max_arrows)
        start_x, step_x, arrow_nx = calculate_sampling(nx, step_size)
        start_y, step_y, arrow_ny = calculate_sampling(ny, step_size)
    else
      start_y, step_y = calculate_start_step(ny, arrow_ny)
      start_y, step_y, arrow_ny = calculate_sampling(ny, step_y)
      start_x, step_x = calculate_start_step(nx, arrow_nx)
      start_x, step_x, arrow_nx = calculate_sampling(nx, step_x)
    end

    I = start_x .+ (0:(arrow_nx - 1)) .* step_x
    J = start_y .+ (0:(arrow_ny - 1)) .* step_y

    Dx = dx * step_x
    Dy = dy * step_y

    if fig === nothing
        fig = Figure(; size=(size_x, size_y), backgroundcolor=:white)
    end

    if colorrange == nothing
        colorrange = :automatic
    end

    if ax == nothing
        ax = Axis(fig[1, 1]; width=size_x, height=size_y)
    end
    hidedecorations!(ax)
    
    heatmap!(ax, xs, ys, mz; alpha=0.5)
    #scatter!(ax, [(x, y) for x in xs for y in ys], color=:white, strokecolor=:black, strokewidth=0.5)

    lengthscale = 0.3 * sqrt(Dx^2 + Dy^2)
    #FIXME: it seems that align=:center does not work well for some situations?
    arrows!(ax, xs[I], ys[J], mx[I, J], my[I, J]; linewidth=2.0, color=:gray36,
            lengthscale=lengthscale, align=:center, colorrange=colorrange)

    return fig
end

"""
    plot_m(sim; k=1, arrows=(-1, -1), figsize=(600, -1))

Create a plotting for given magnetization. 
  `sim` should be a Sim Object.
  `k` indicates the layer index (starting from 1) 
  `arrows` is the number of arrows, should be a Tuple of integers. By default, arrows=(-1, -1).
  `figsize` should be a Tuple of integers, for example, figsize=(500, 400) or figsize=(500, -1).
"""
function JuMag.plot_m(sim::JuMag.AbstractSim; k=1, arrows=(-1, -1), figsize=(600, -1))
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    m = Array(sim.spin)
    m = reshape(m, 3, nx, ny, nz)
    fig = plot_m(m; k=k, dx=mesh.dx, dy=mesh.dy, arrows=arrows, figsize=figsize)
    return fig
end

"""
  ovf2png(ovf_name, output=nothing; k=1, arrows=(-1, -1), figsize=(500, -1))

Create a png from the given ovf file.
`k` indicates the layer index (starting from 1) 
`arrows` is the number of arrows, should be a Tuple of integers. By default, arrows=(-1, -1).
`figsize` should be a Tuple of integers, for example, figsize=(500, 400) or figsize=(500, -1).
"""
function JuMag.ovf2png(ovf_name, output=nothing; k=1, arrows=(-1, -1), figsize=(500, -1))
    if output === nothing
        output = endswith(ovf_name, ".ovf") ? ovf_name[1:(end - 4)] : ovf_name
    end
    ovf = read_ovf(ovf_name)
    spin = reshape(ovf.data, 3, ovf.xnodes, ovf.ynodes, ovf.znodes)
    fig = JuMag.plot_m(spin; dx=ovf.xstepsize, dy=ovf.ystepsize, k=k, arrows=arrows,
                 figsize=figsize)
    save(output * ".png", fig)
    return fig
end

"""
  jdl2movie(jdl_file; k = 1, component='z', framerate=30, rm_png=false, output=nothing)

Create a moive from the given jdl2 file. `k` indicates the layer index (starting from 1), 
`component` refers to the used component ('x', 'y' or 'z') for plotting, `r` is frame rate, 
output is the filename of the video and the support formats are 'mp4', 'avi' and 'gif'.
"""
function JuMag.jdl2movie(jdl_file; k = 1, component='z', figsize=(500, -1), framerate=12, output=nothing)
  if output===nothing
    base_name = jdl_file[1:length(jdl_file)-5]
    output = @sprintf("%s.mp4", base_name)
  end

  data = JLD2.load(jdl_file)
  steps = data["steps"]
  save_m_every = data["save_m_every"]
  nx, ny, nz = data["mesh/nx"], data["mesh/ny"], data["mesh/nz"]
  (k < 1) && (k = 1)
  (k > nz) && (k = nz)
  m_index = component - 'x' + 1
  if save_m_every < 0
    @info @sprintf("save_m_every is %d, which is negative, exiting~", save_m_every)
    return
  end

  dx, dy, dz = data["mesh/dx"], data["mesh/dy"], data["mesh/dz"]

  size_x = figsize[1]
  size_y = figsize[2]
  if (size_y < 0)
      aspect_ratio = ny * dy / (nx * dx)
      size_y = Int(ceil(size_x * aspect_ratio))
  end

  fig = Figure(; size=(size_x, size_y), backgroundcolor=:white)

  ax = Axis(fig[1, 1]; width=size_x, height=size_y)
  hidedecorations!(ax)

  function update_function(i)
    index = @sprintf("m/%d", i)
    m = reshape(data[index], 3, nx, ny, nz)
    plot_m(m, dx=dx, dy=dy, k=k, fig=fig, ax=ax, colorrange=colorrange=[-1, 1])
  end

  record(update_function, fig, output, 1:save_m_every:steps; framerate = framerate)
end

end
