using PyCall

function m2rgb(mx,my,mz,rotation)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    plt = pyimport("matplotlib.pyplot")
    colorsys = pyimport("colorsys")

    nx,ny = size(mx)
    data_hsv = zeros(nx,ny,3)
    for i =1:nx,j=1:ny
        data_hsv[i,j,1] = atan(my[i,j],mx[i,j]) + pi + rotation
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

function plot_m(fname, mx, my, mz, component, rotation, quiver_args)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")

    if component == "mz"
        plt.imshow(np.transpose(mz), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "mx"
        plt.imshow(np.transpose(mx), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "my"
        plt.imshow(np.transpose(my), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "all"
        max_value = maximum([maximum(mx),maximum(my),maximum(mz)])
        mx /= max_value
        my /= max_value
        mz /= max_value
        data_rgb = m2rgb(mx,my,mz,rotation)
        ls = mcolors.LightSource()
        plt.imshow(data_rgb, origin="lower", interpolation="bicubic")
    end

    plt.xticks([])
    plt.yticks([])

    if quiver_args[:quiver]
        x,y = size(mx)
        interval = quiver_args[:interval]
        X,Y = np.meshgrid([1:interval:x],[1:interval:y])
        norm_constant = max(maximum(abs.(mx)), maximum(abs.(my)))
        mx ./= norm_constant
        my ./= norm_constant
        plt.axes().set_aspect("equal")
        prog = "plt.quiver(X,Y,transpose(mx[1:interval:end,1:interval:end]),transpose(my[1:interval:end,1:interval:end])"
        delete!(quiver_args, :quiver)
        delete!(quiver_args, :interval)
        for (key,value) in quiver_args
            if typeof(value) == String
                prog *= ", $key = \"$value\""
            else 
                prog *= ", $key = $value" 
            end
        end

        prog *= ")"
        print(prog)
        eval(Meta.parse(prog))
    end

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=dirname(fname)

    if !isdir(path)
        mkpath(path)
    end
    
    plt.tight_layout()
    plt.savefig(joinpath(path,basename(fname))*"_slice_$component.png",dpi=300)
    plt.close()
end


"""
  plotOVF(fname; slice="center", component = "mz", cmap="coolwarm", rotation = 0, quiver=false, quiver_interval=5, quiver_size=30, quiver_color="w")

Plot a certain slice of the given ovf file.

Parameters
-----------------
fname: file name
slice: layer number.
component: Can be chosen from "mx","my","mz" and "all".
rotation: Rotation angle of the hsv color space in degree.
quiver_args: A dictionary used in python plt.quiver function.
    For example: quiver_args=Dict(:quiver => true, :interval =>1, color => "w", scale =>20) will execute 
                plt.quiver(mx[::interval,::interval],my[::interval,::interval],color="w", scale=20)
"""
function plotOVF(fname; slice="center", component = "mz", cmap="coolwarm", rotation = 0, 
    quiver_args::Dict=Dict(:quiver =>false))
    ovf = read_ovf(fname)
    m = ovf.data
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    b = reshape(m, (3,nx,ny,nz))
    if slice == "center"
        slice = ceil(Int, nz/2)
    end
    mx, my, mz = b[1, :, :, slice], b[2, :, :, slice], b[3, :, :, slice]

    #= quiver_args = Dict(:quiver => quiver, :interval => quiver_interval,
                        :scale => quiver_size, :color => quiver_color) =#

    plot_m(fname, mx, my, mz, component, rotation, quiver_args)
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
