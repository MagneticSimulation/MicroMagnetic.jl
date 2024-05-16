
module CairoMakieExt

using JuMag
using CairoMakie

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
    plot_m(spin; dx=1.0, dy=1.0, k=1, arrows=(-1, 1), figsize=500)

Create a plotting for given magnetization. 

  `spin` should be an array with dimension (3, nx, ny, nz)
  `k` indicates the layer index (starting from 1) 
  `arrows` is the number of arrows, should be a Tuple of integers. By default, arrows=(-1, -1).
  `figsize` should be a Tuple of integers, for example, figsize=(500, 400) or figsize=(500, -1).
"""
function JuMag.plot_m(spin; dx=1.0, dy=1.0, k=1, arrows=(-1, -1), figsize=(500, -1))
    (_, nx, ny, nz) = size(spin)
    scale_factor = 10^floor(log10(dx))
    dx = dx / scale_factor
    dy = dy / scale_factor
    xs = [i * dx for i in 1:nx]
    ys = [j * dy for j in 1:ny]

    mx = spin[1, :, :, k]
    my = spin[2, :, :, k]
    mz = spin[3, :, :, k]

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
      start_x, step_x = calculate_start_step(nx, arrow_nx)
    end

    I = start_x .+ (0:(arrow_nx - 1)) .* step_x
    J = start_y .+ (0:(arrow_ny - 1)) .* step_y

    Dx = dx * step_x
    Dy = dy * step_y

    fig = Figure(; size=(size_x, size_y), backgroundcolor=:white)

    ax = Axis(fig[1, 1]; width=size_x, height=size_y)
    hidedecorations!(ax)

    heatmap!(ax, xs, ys, mz; alpha=0.5)
    #scatter!(ax, [(x, y) for x in xs for y in ys], color=:white, strokecolor=:black, strokewidth=0.5)

    lengthscale = 0.3 * sqrt(Dx^2 + Dy^2)
    #FIXME: it seems that align=:center does not work well for some situations?
    arrows!(ax, xs[I], ys[J], mx[I, J], my[I, J]; linewidth=2.0, color=:gray36,
            lengthscale=lengthscale, align=:center)

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

end
