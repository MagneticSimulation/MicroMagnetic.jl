using PyCall
using MicroMag


"""
    mag_to_rgb(mx,my,mz,rotation)

Compute rgb array  

Input:  m : Tuple(mx, my, mz), where mx is a 2d array sized (nx,ny)
        ra : rgb rotation angle in rad

output: rgb array with size (3, nx, ny)

"""

function mag_to_rgb(m::Tuple, ra::Number)
    colorsys = pyimport("colorsys")

    (mx, my, mz) = m

    nx,ny = size(mx)
    data_hsv = zeros(nx,ny,3)
    for i =1:nx,j=1:ny
        data_hsv[i,j,1] = atan(my[i,j],mx[i,j]) + pi + ra
        data_hsv[i,j,1] = mod(data_hsv[i,j,1], 2*pi)

        data_hsv[i,j,2] = sqrt(mx[i,j]^2+my[i,j]^2)
        data_hsv[i,j,3] = mz[i,j]*0.5 + 0.5
    end
    data_hsv[:,:,1] /= 2*pi
    data_rgb = zeros(ny,nx,3)
    for i =1:nx,j=1:ny
        data_rgb[j,i,:] .= colorsys.hsv_to_rgb(data_hsv[i,j,1],data_hsv[i,j,2],data_hsv[i,j,3])
    end
    return data_rgb
end

function show_mag(m::Union{Tuple, Array{T, 2}}, fname::String; ra=0) where T <: AbstractFloat
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")

    if m isa Array
        plt.imshow(np.transpose(m), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif m isa Tuple
        max_value = maximum([maximum(m[1]), maximum(m[2]), maximum(m[3])])
        mx = m[1]/max_value
        my = m[2]/max_value
        mz = m[3]/max_value
        data_rgb = mag_to_rgb((mx,my,mz), ra)
        ls = mcolors.LightSource()
        plt.imshow(data_rgb, origin="lower", interpolation="bicubic")
    else
        @error("Input must be a 2d array or a tuple of three 2d array!")
    end

    plt.xticks([])
    plt.yticks([])

    if endswith(fname, ".ovf")
        fname = fname[1:end-4]
    end

    path=dirname(fname)

    if !isdir(path)
        mkpath(path)
    end

    plt.tight_layout()
    plt.savefig(joinpath(path,basename(fname))*".png",dpi=300)
    plt.close()
end

"""
    show_mag(fname::String; layer=0, rgb=false, component = 3, ra = 0)

Show the magnetization of an ovf file.

Parameters:

    fname : .ovf file name including relative path.

Optional:

    layer : Which layer to be shown. Integer from 1 to nz.

    component : 1 or 2 or 3, representing mx, my, and mz respectively.

    rgb : If true, this function will use the rgb colorwheel.

    ra : rgb rotation angle in rad.

Example:
```julia
    show_mag("example.ovf", rgb=true)
```
"""
function show_mag(fname::String; layer=0, rgb=false, component = 3, ra = 0)

    ovf = read_ovf(fname)
    m = ovf.data
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    b = reshape(m, (3,nx,ny,nz))
    if layer == 0
        layer = ceil(Int, nz/2)
    end
    mx, my, mz = b[1, :, :, layer], b[2, :, :, layer], b[3, :, :, layer]
    
    if rgb
        show_mag((mx,my,mz), fname, ra=ra)
    else
        show_mag(b[component, :, :, layer], fname)
    end
end

"""
  plot_ovf_projection(fname; tilt_angle=0, tilt_axis="alpha", component = "mz", cmap="coolwarm", rotation = 0, quiver=false, quiver_interval=5, quiver_size=30, quiver_color="w", N=128)

Plot the averaged magnetization from the given ovf file.

Based on radon transform.

Parameters
--------------------
tilt_axis: Can be "alpha"(x-axis) or "beta"(y-axis).
tilt_angle: Tilt angle in degree.

Other parameters are same with function "plot_ovf_slice"
"""

function plot_ovf_projection(fname; tilt_angle=0, tilt_axis="alpha", component = "mz", cmap="coolwarm", rotation = 0, 
    quiver_args::Dict=[:quiver =>false], N=128)

    ovf = read_ovf(fname)
    m = reshape(ovf.data, (3, ovf.xnodes, ovf.ynodes, ovf.znodes))
    m = vector_padding(m, N,N,N)
    mp = vector_field_projection(m, tilt_angle, tilt_axis)
    mx, my, mz = mp[1,:,:], mp[2,:,:], mp[3,:,:]

    plot_m(fname, mx, my, mz, component, rotation, quiver_args)   
end
