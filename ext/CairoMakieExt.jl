
module CairoMakieExt

using Printf
using CairoMakie
using MicroMagnetic

function calculate_sampling_parameters(nx::Int, ny::Int, arrows::Tuple{Int,Int})
    """
    Calculate sampling parameters with fractional start positions
    Ensures symmetry and uniform spacing
    """
    arrow_nx, arrow_ny = arrows
    
    function symmetric_sampling(n_total::Int, n_samples::Int)
        # Single point case - place at center
        if n_samples <= 1
            return (n_total + 1) / 2, 1.0, 1
        end

        step = n_total / n_samples
        start = step / 2 + 0.5
        return start, step, n_samples
    end
    
    function adaptive_sampling(n_total::Int, reference_step::Real)        
        # Calculate number of samples
        n_samples = max(1, round(Int, n_total / reference_step))
        step = n_total / n_samples
        start = 0.5 + step/2
        return start, step, n_samples
    end
    
    # Different sampling modes based on arrow configuration
    if arrows[1] < 0 && arrows[2] > 0
        # Mode 1: Fixed y samples, adaptive x
        start_y, step_y, arrow_ny = symmetric_sampling(ny, arrow_ny)
        start_x, step_x, arrow_nx = adaptive_sampling(nx, step_y)
        
    elseif arrows[1] > 0 && arrows[2] < 0
        # Mode 2: Fixed x samples, adaptive y
        start_x, step_x, arrow_nx = symmetric_sampling(nx, arrow_nx)
        start_y, step_y, arrow_ny = adaptive_sampling(ny, step_x)
        
    elseif arrows[1] < 0 && arrows[2] < 0
        # Mode 3: Both directions adaptive, same density
        max_arrows = 40
        _, step_size, _ = symmetric_sampling(max(nx, ny), max_arrows)
        start_x, step_x, arrow_nx = adaptive_sampling(nx, step_size)
        start_y, step_y, arrow_ny = adaptive_sampling(ny, step_size)
        
    else
        # Mode 4: Both directions fixed
        start_x, step_x, arrow_nx = symmetric_sampling(nx, arrow_nx)
        start_y, step_y, arrow_ny = symmetric_sampling(ny, arrow_ny)
    end
    
    return start_x, step_x, arrow_nx, start_y, step_y, arrow_ny
end

function bilinear_interpolation_grid(data::AbstractMatrix, x_coords, y_coords)
    nx, ny = size(data, 1), size(data, 2)
    new_nx = length(x_coords)
    new_ny = length(y_coords)
    result = zeros(new_nx, new_ny)
    
    for j in 1:new_nx
        x = x_coords[j]
        ix = clamp(floor(Int, x), 1, nx - 1)
        tx = x - ix
        
        for i in 1:new_ny
            y = y_coords[i]
            iy = clamp(floor(Int, y), 1, ny - 1)
            ty = y - iy

            bottom_val = data[ix, iy] * (1 - tx) + data[ix+1, iy] * tx
            top_val = data[ix, iy+1] * (1 - tx) + data[ix+1, iy+1] * tx
            result[j, i] = bottom_val * (1 - ty) + top_val * ty
        end
    end
    
    return result
end

"""
    MicroMagnetic.plot_m(spin; dx=1.0, dy=1.0, k=1, component='z', arrows=(-1, -1), figsize=(500, -1), fig=nothing, ax=nothing, kwargs...)

Create a plot for the given magnetization.

# Arguments
- `spin::Array{Float64, 4}`: An array with dimensions (3, nx, ny, nz) representing the magnetization.

# Keyword Arguments
- `dx::Float64`: The spacing in the x-direction (default: 1.0).
- `dy::Float64`: The spacing in the y-direction (default: 1.0).
- `k::Int`: The layer index to plot (starting from 1) (default: 1).
- `component::Char`: The magnetization component to plot ('x', 'y', or 'z') (default: 'z').
- `arrows::Tuple{Int, Int}`: The number of arrows to plot, specified as a tuple. By default, arrows=(-1, -1), which auto-scales the number of arrows. If arrows is set to `nothing`, no arrows will be shown.
- `figsize::Tuple{Int, Int}`: The size of the figure, specified as a tuple (width, height). For example, figsize=(500, 400) or figsize=(500, -1) where -1 means auto-scaled height (default: (500, -1)).
- `fig`: An existing figure to plot on. If nothing, a new figure is created (default: nothing).
- `ax`: An existing axis to plot on. If nothing, a new axis is created (default: nothing).
- `kwargs...`: Additional keyword arguments that are passed to `heatmap!`.

# Supported `heatmap!` Arguments
- `alpha`: Transparency level of the heatmap.
- `colormap`: A list or function that defines the color mapping for the heatmap [Colors](https://docs.makie.org/v0.21/explanations/colors).
- `colorrange`: The color range of the heatmap.
- `interpolate`: Whether to interpolate between values.
- `colorscale`: A scale factor for the color.
- `transparency`: Whether to apply transparency to the heatmap.

For more details on the supported arguments, please refer to the Makie documentation: [Makie Heatmap Reference](https://docs.makie.org/v0.21/reference/plots/heatmap).

# Returns
- `fig`: The figure containing the plot.

# Examples
```julia
spin = randn(3, 10, 10, 5)  # Example spin data
# Creates a plot with default settings
MicroMagnetic.plot_m(spin, colorrange=[-1, 1]) 

# Creates a plot for the x-component of the second layer with custom settings
MicroMagnetic.plot_m(spin, dx=0.5, dy=0.5, k=2, component='x', arrows=(5, 5), figsize=(600, 400))
```
"""
function MicroMagnetic.plot_m(spin; dx=1.0, dy=1.0, k=1, component='z', arrows=(-1, -1),
                              figsize=(500, -1), fig=nothing, ax=nothing, kwargs...)
    (_, nx, ny, nz) = size(spin)
    scale_factor = 10^floor(log10(dx))
    dx = dx / scale_factor
    dy = dy / scale_factor
    xs = [i * dx for i in 1:nx]
    ys = [j * dy for j in 1:ny]

    mx = spin[1, :, :, k]
    my = spin[2, :, :, k]
    mz = spin[3, :, :, k]
    lml = sqrt.(mx .^ 2 .+ my .^ 2 .+ mz .^ 2)
    mx[lml .< 0.1] .= NaN
    my[lml .< 0.1] .= NaN
    mz[lml .< 0.1] .= NaN

    size_x = figsize[1]
    size_y = figsize[2]
    if (size_y < 0)
        aspect_ratio = ny * dy / (nx * dx)
        size_y = Int(ceil(size_x * aspect_ratio))
    end

    show_arrow = true
    if arrows === nothing
        arrows = (-1, -1)
        show_arrow = false
    end

    start_x, step_x, arrow_nx, start_y, step_y, arrow_ny = calculate_sampling_parameters(nx, ny, arrows)

    I = start_x .+ (0:(arrow_nx - 1)) .* step_x
    J = start_y .+ (0:(arrow_ny - 1)) .* step_y

    Dx = dx * step_x
    Dy = dy * step_y

    if fig === nothing
        fig = Figure(; size=(size_x, size_y), backgroundcolor=:white)
    end

    if ax === nothing
        ax = Axis(fig[1, 1]; width=size_x, height=size_y)
    end
    hidedecorations!(ax)

    if component == 'x'
        mm = mx
    elseif component == 'y'
        mm = my
    else
        mm = mz
    end

    valid_heatmap_args = [:alpha, :colormap, :colorrange, :interpolate, :colorscale,
                          :transparency]
    filtered_kwargs = Dict(kw => val for (kw, val) in kwargs if kw in valid_heatmap_args)

    heatmap!(ax, xs, ys, mm; filtered_kwargs...)

    if show_arrow
        lengthscale = 0.5 * sqrt(Dx^2 + Dy^2)
        mx_interp = bilinear_interpolation_grid(mx, I, J)
        my_interp = bilinear_interpolation_grid(my, I, J)
        arrows2d!(ax, I.*dx, J.*dy, mx_interp, my_interp; tipwidth=11, shaftwidth=5, color=:gray36,
               lengthscale=lengthscale, align=:center)
    end
    return fig
end

"""
    MicroMagnetic.plot_m(sim::MicroMagnetic.AbstractSim; kwargs...)

Create a plot for the given magnetization in a simulation object.

# Arguments
- `sim::MicroMagnetic.AbstractSim`: A simulation object containing the magnetization data and mesh information.

# Keyword Arguments
This function forwards all keyword arguments to `MicroMagnetic.plot_m`. Refer to `MicroMagnetic.plot_m` for detailed descriptions of the keyword arguments.

# Returns
- `fig`: The figure containing the plot.

# Examples
```julia
sim = MicroMagnetic.create_simulation()  # Example simulation object
MicroMagnetic.plot_m(sim)
# Creates a plot with default settings

MicroMagnetic.plot_m(sim, k=2, component='x', arrows=(5, 5), figsize=(600, 400))
# Creates a plot for the x-component of the second layer with custom settings
```
"""
function MicroMagnetic.plot_m(sim::MicroMagnetic.AbstractSim; kwargs...)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    m = Array(sim.spin)
    m = reshape(m, 3, nx, ny, nz)
    fig = MicroMagnetic.plot_m(m; dx=mesh.dx, dy=mesh.dy, kwargs...)
    return fig
end

"""
    ovf2png(ovf_name, output=nothing; k=1, arrows=(-1, -1), figsize=(500, -1))

Create a png from the given ovf file.

- `k` indicates the layer index (starting from 1) 
- `arrows` is the number of arrows, should be a Tuple of integers. By default, arrows=(-1, -1).
- `figsize` should be a Tuple of integers, for example, figsize=(500, 400) or figsize=(500, -1).
"""
function MicroMagnetic.ovf2png(ovf_name, output=nothing; k=1, arrows=(-1, -1),
                               figsize=(500, -1))
    if output === nothing
        output = endswith(ovf_name, ".ovf") ? ovf_name[1:(end - 4)] : ovf_name
    end
    ovf = read_ovf(ovf_name)
    spin = reshape(ovf.data, 3, ovf.xnodes, ovf.ynodes, ovf.znodes)
    fig = MicroMagnetic.plot_m(spin; dx=ovf.xstepsize, dy=ovf.ystepsize, k=k, arrows=arrows,
                               figsize=figsize)
    save(output * ".png", fig)
    return fig
end

"""
   ovf2movie(folder; framerate=12, output=nothing, kwargs...)

Create a moive from the given folder, 

# Keyword Arguments
This function forwards all keyword arguments to `MicroMagnetic.plot_m`. Refer to `MicroMagnetic.plot_m` for detailed descriptions of the keyword arguments.

`output` is the filename of the video and the support formats are 'mp4', 'avi' and 'gif'.

# Examples
```julia
ovf2movie("std4_LLG"; output="std4.gif", component='x');
```
"""
function MicroMagnetic.ovf2movie(folder; framerate=12, output=nothing, figsize=(500, -1),
                                 kwargs...)

    if output === nothing
        output = @sprintf("%s.gif", basename(path))
    end

    files = readdir(folder)
    ovf_files = filter(f -> startswith(f, "m_") && endswith(f, ".ovf"), files)
    sorted_files = sort(ovf_files, by = f -> parse(Int, match(r"m_(\d+)\.ovf", f).captures[1]))
    
    ovf = read_ovf(folder*"/"*sorted_files[1])
    dx, dy = ovf.xstepsize, ovf.ystepsize
    nx, ny = ovf.xnodes, ovf.ynodes
    
    size_x = figsize[1]
    size_y = figsize[2]
    if (size_y < 0)
        aspect_ratio = ny * dy / (nx * dx)
        size_y = Int(ceil(size_x * aspect_ratio))
    end

    fig = Figure(; size=(size_x, size_y), backgroundcolor=:white)

    ax = Axis(fig[1, 1]; width=size_x, height=size_y)
    hidedecorations!(ax)

    function update_function(ovf_name)
        ovf = read_ovf(folder*"/"*ovf_name)
        spin = reshape(ovf.data, 3, ovf.xnodes, ovf.ynodes, ovf.znodes)
        empty!(ax)
        return MicroMagnetic.plot_m(spin; dx=ovf.xstepsize, dy=ovf.ystepsize, fig=fig, ax=ax, kwargs...)
    end

    record(update_function, fig, output, sorted_files; framerate=framerate)
    
    return output
end


"""
    plot_ts(filename::String, keys::Vector{String}; 
            x_key::String="time", x_unit::Real=1e9, y_units::Vector=[],
            size=(400, 280), title="", xlabel="", ylabel="", 
            plot_type=:scatterlines, markersize=6, legend_position=:rt,
            colors=nothing, linestyles=nothing, transparency=false)

General time series plotting function

# Arguments
- `filename`: Data file name
- `keys`: Data column keys to plot
- `x_key`: Time column key, default is "time"
- `x_unit`: Time unit conversion factor, default is 1e9 (seconds to nanoseconds)
- `y_units`: Y-axis data unit conversion factors for each key, same order as `keys`
- `plot_type`: Plot type, `:lines`, `:scatter`, `:scatterlines`, `:linesmarkers`
- `legend_position`: Legend position, `:rt`(top right), `:lt`(top left), `:rb`(bottom right), `:lb`(bottom left)
- `colors`: Custom color sequence
- `linestyles`: Custom line style sequence
- `transparency`: Whether to use transparent background

# Examples
```julia
# Plot magnetization components
plot_ts("std4_llg.txt", ["m_x", "m_y", "m_z"]; 
        xlabel="Time (ns)", ylabel="m", plot_type=:scatter)

# Plot vortex center with unit conversion
plot_ts("std5_llg.txt", ["cx", "cy"]; 
        xlabel="Time (ns)", ylabel="Vortex center (nm)", 
        y_units=[1e9, 1e9], plot_type=:scatterlines)

# Custom colors and line styles
plot_ts("data.txt", ["A", "B", "C"]; 
        colors=[:red, :blue, :green],
        linestyles=[:solid, :dash, :dot])
```
"""
function MicroMagnetic.plot_ts(filename::String, keys::Vector{String}; 
                x_key::String="time", x_unit::Real=1e9, y_units::Vector=[],
                size=(400, 280), title="", xlabel="", ylabel="", 
                plot_type=:scatterlines, markersize=6, legend_position=:rt,
                colors=nothing, linestyles=nothing, transparency=false)
    
    # Load data
    data, unit = MicroMagnetic.read_table(filename)
    
    # Create figure
    bg_color = transparency ? :transparent : :white
    fig = Figure(; size=size, backgroundcolor=bg_color)
    ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=title, backgroundcolor=bg_color)
    
    # Time data
    time_data = data[x_key] .* x_unit
    
    # Default colors and line styles
    default_colors = [:slateblue1, :sienna1, :seagreen2, :tomato1, :gold1,
                      :darkorchid1, :deepskyblue2, :coral1, :limegreen, :hotpink,
                      :royalblue1, :darkorange1, :mediumseagreen, :violetred1, :steelblue2,
                      :chocolate1, :mediumpurple1, :darkturquoise, :orangered, :mediumvioletred]
    default_linestyles = [:solid, :dash, :dot, :dashdot]
    
    colors = isnothing(colors) ? default_colors : colors
    linestyles = isnothing(linestyles) ? default_linestyles : linestyles
    
    # Plot each data series
    for (i, key) in enumerate(keys)
        # Get data and apply unit conversion
        y_data = data[key]
        if !isempty(y_units) && i <= length(y_units)
            y_data = y_data .* y_units[i]
        end
        
        # Select color and line style
        color = colors[mod1(i, length(colors))]
        linestyle = linestyles[mod1(i, length(linestyles))]
        
        # Plot according to plot type
        if plot_type == :lines
            lines!(ax, time_data, y_data; color=color, linestyle=linestyle, label=key)
        elseif plot_type == :scatter
            scatter!(ax, time_data, y_data; color=color, markersize=markersize, label=key)
        elseif plot_type == :scatterlines
            scatterlines!(ax, time_data, y_data; color=color, markersize=markersize, 
                         linestyle=linestyle, label=key)
        elseif plot_type == :linesmarkers
            lines!(ax, time_data, y_data; color=color, linestyle=linestyle, label=key)
            scatter!(ax, time_data, y_data; color=color, markersize=markersize)
        else
            error("Unknown plot_type: $plot_type. Use :lines, :scatter, :scatterlines, or :linesmarkers")
        end
    end
    
    # Add legend
    if length(keys) > 1
        legend_pos = if legend_position == :rt
            :rt
        elseif legend_position == :lt
            :lt
        elseif legend_position == :rb
            :rb
        elseif legend_position == :lb
            :lb
        else
            :rt
        end
        axislegend(ax; position=legend_pos)
    end
    
    return fig
end

function MicroMagnetic.plot_voronoi(grain_ids, points; dx=2, dy=2, output="voronoi.png")
    
    nx, ny = size(grain_ids)
    Lx, Ly = nx*dx, ny*dy

    fig = Figure()
    ax = Axis(fig[1, 1])
    heatmap!(ax, 0:Lx, 0:Ly, grain_ids, colormap=Reverse(:viridis))
    
    for i in 1:length(points)
        p = points[i]
        text!(ax,  @sprintf("%d", i), 
            position = (round(Int, p[1]), round(Int, p[2])), 
            align = (:center, :center),
            fontsize = 14
        )
    end
    
    save(output, fig)
end

end
