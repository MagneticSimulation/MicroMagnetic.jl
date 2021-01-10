#using PyCall

function M2RGB(mx,my,mz)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")
    colorsys = pyimport("colorsys")

    nx,ny = size(mx)
    data_hsv = zeros(nx,ny,3)
    for i =1:nx,j=1:ny
        data_hsv[i,j,1] = atan(my[i,j],mx[i,j])
        data_hsv[i,j,1] = data_hsv[i,j,1]/(2*pi) + 0.25
        if data_hsv[i,j,1]<0
            data_hsv[i,j,1] = data_hsv[i,j,1]+1
        end

        data_hsv[i,j,2] = sqrt(mx[i,j]^2+my[i,j]^2)
        data_hsv[i,j,3] = -1*mz[i,j]^2+1
    end
    data_rgb = zeros(nx,ny,3)
    for i =1:nx,j=1:ny
        data_rgb[i,j,:] .= colorsys.hsv_to_rgb(data_hsv[i,j,1],data_hsv[i,j,2],data_hsv[i,j,3])
    end
    return data_rgb
end
"""
    plotOVF(fname;axis=ez,slice=-1,style = "mz",norm=nothing, cmap="coolwarm")

plot the projection magnetization from a .ovf file.

axis: projection axis

slice: Chosen layer. Set -1 to get projection from all layers by default.

style: Chosen component along axis:

    "mx", "my", "mz": single projection component.

    "rgb": hsv color space

the normalization can be set by using norm=[vmin,vmax]
"""
function plotOVF(fname;axis=ez,slice=-1,style = "mz",norm=nothing, cmap="coolwarm",quiver=false,quiver_interval=3,quiver_size=45)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    plt = pyimport("matplotlib.pyplot")
    mcolors = pyimport("matplotlib.colors")
    colorsys = pyimport("colorsys")

    ovf=read_ovf(fname)
    mxp,myp,mzp = ovf_projection(ovf,axis=axis,slice=slice)
    if style == "mz"
        plt.imshow(np.transpose(mzp), origin="lower", cmap="coolwarm")
        #plt.colorbar(orientation="horizontal",norm = norm)
        plt.colorbar()
    elseif style == "mx"
        plt.imshow(np.transpose(mxp), origin="lower", cmap="coolwarm")
         plt.colorbar(norm = norm)
    elseif style == "my"
        plt.imshow(np.transpose(myp), origin="lower", cmap="coolwarm")
         plt.colorbar(norm = norm)
    elseif style == "mxy"
        plt.imshow(np.transpose(sqrt.(mxp.^2+myp.^2)), origin="lower", cmap="coolwarm")
        plt.colorbar(norm = norm)
    elseif style == "rgb"
        nz = ovf.znodes
        mx,my,mz = mxp./nz,myp./nz,mzp./nz
        data_rgb = M2RGB(mx,my,mz)
        ls = mcolors.LightSource()
        a,b,c = size(data_rgb)
        data_rgb1 = zeros(b,a,c)
        for i =1:3
            data_rgb1[:,:,i] .= np.transpose(data_rgb[:,:,i])
        end
        plt.imshow(data_rgb1, origin="lower", interpolation="bicubic")
    end

    if norm != nothing
        #norm = mpl.colors.Normalize(vmin=norm[1],vmax=norm[2])
        plt.clim(norm[1],norm[2])
    end

    plt.axis("off")

    if quiver
        x,y = size(mxp)
        X,Y = np.meshgrid([1:quiver_interval:x],[1:quiver_interval:y])
        plt.axes().set_aspect("equal")
        plt.quiver(X,Y,np.transpose(mxp[1:quiver_interval:end,1:quiver_interval:end]),np.transpose(myp[1:quiver_interval:end,1:quiver_interval:end]),
            pivot="mid",scale=quiver_size)
    end

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"m")
    if !isdir(path)
        mkpath(path)
    end
    plt.savefig(joinpath(path,basename(fname))*"_$style.jpg",dpi=300,bbox_inches="tight",pad_inches=0.0)
    plt.close()
end

function plot_folder(folder;axis=ez,slice=-1,style = "mz",norm=nothing, cmap="coolwarm",quiver=false,quiver_interval=3,quiver_size=45)
    files = readdir(folder)
    for f in files
        if endswith(f,".ovf")
            plotOVF(joinpath(folder,f),
                axis=axis,slice=slice,style = style,norm=norm, cmap=cmap,quiver=quiver,quiver_interval=quiver_interval,quiver_size=quiver_size)
        end
    end
end
