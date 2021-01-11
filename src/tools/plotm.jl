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

function plot_ovf_slice(fname; slice="center", component = "mz", cmap="coolwarm", rotation = 0, quiver=false, quiver_interval=5, quiver_size=30, quiver_color="w")
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")

    ovf = read_ovf(fname)
    m = ovf.data
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    b = reshape(m, (3,nx,ny,nz))
    if slice == "center"
        slice = round(Int, nz/2)
    end
    mx, my, mz = b[1, :, :, slice], b[2, :, :, slice], b[3, :, :, slice]

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
        biggest = maximum([maximum(mx),maximum(my),maximum(mz)])
        mx /= biggest
        my /= biggest
        mz /= biggest
        data_rgb = m2rgb(mx,my,mz,rotation)
        ls = mcolors.LightSource()
        plt.imshow(data_rgb, origin="lower", interpolation="bicubic")
    end

    plt.xticks([])
    plt.yticks([])

    if quiver
        x,y = size(mx)
        X,Y = np.meshgrid([1:quiver_interval:x],[1:quiver_interval:y])
        plt.axes().set_aspect("equal")
        plt.quiver(X,Y,np.transpose(mx[1:quiver_interval:end,1:quiver_interval:end]),np.transpose(my[1:quiver_interval:end,1:quiver_interval:end]),
            pivot="mid",scale=quiver_size,color=quiver_color)
    end

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"m")
    if !isdir(path)
        mkpath(path)
    end
    plt.tight_layout()
    plt.savefig(joinpath(path,basename(fname))*"_$component.png",dpi=300)
    plt.close()

end

function plot_ovf_projection(fname; angle=0, tilt_axis="alpha", component = "mz", cmap="coolwarm", rotation = 0, quiver=false, quiver_interval=5, quiver_size=30, quiver_color="w", N=128)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")

    ovf = read_ovf(fname)
    mxp,myp,mzp = radon_transform_ovf(ovf, angle, tilt_axis, N=N)
    if component == "mz"
        plt.imshow(np.transpose(mzp), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "mx"
        plt.imshow(np.transpose(mxp), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "my"
        plt.imshow(np.transpose(myp), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "mxy"
        plt.imshow(np.transpose(sqrt.(mxp.^2+myp.^2)), origin="lower", cmap="coolwarm")
        plt.colorbar()
    elseif component == "all"
        biggest = maximum([maximum(mxp),maximum(myp),maximum(mzp)])
        mxp /= biggest
        myp /= biggest
        mzp /= biggest
        data_rgb = m2rgb(mxp,myp,mzp,rotation)
        ls = mcolors.LightSource()
        plt.imshow(data_rgb, origin="lower", interpolation="bicubic")
    end

    plt.xticks([])
    plt.yticks([])

    if quiver
        x,y = size(mxp)
        X,Y = np.meshgrid([1:quiver_interval:x],[1:quiver_interval:y])
        plt.axes().set_aspect("equal")
        plt.quiver(X,Y,np.transpose(mxp[1:quiver_interval:end,1:quiver_interval:end]),np.transpose(myp[1:quiver_interval:end,1:quiver_interval:end]),
            pivot="mid",scale=quiver_size,color=quiver_color)
    end

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"m")
    if !isdir(path)
        mkpath(path)
    end
    plt.tight_layout()
    plt.savefig(joinpath(path,basename(fname))*"_$component.png",dpi=300)
    plt.close()
end
