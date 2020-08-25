using PyCall

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
        data_hsv[i,j,1] = data_hsv[i,j,1]/(2*pi) + 0.5

        data_hsv[i,j,2] = sqrt(mx[i,j]^2+my[i,j]^2)
        data_hsv[i,j,3] = mz[i,j]*0.5 + 0.5
    end
    data_rgb = zeros(ny,nx,3)
    for i =1:nx,j=1:ny
        data_rgb[j,i,:] .= colorsys.hsv_to_rgb(data_hsv[i,j,1],data_hsv[i,j,2],data_hsv[i,j,3])
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
function plotOVF(fname;axis=ez,slice=-1,style = "mz",norm=nothing, cmap="coolwarm")
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    plt = pyimport("matplotlib.pyplot")

    if norm != nothing
        norm = mpl.colors.Normalize(vmin=norm[1],vmax=norm[2])
    end

    mxp,myp,mzp = ovf_projection(fname,axis=axis,slice=slice)
    if style == "mz"
        plt.imshow(np.transpose(mzp), origin="lower", cmap="coolwarm")
        #plt.colorbar(orientation="horizontal",norm = norm)
        plt.colorbar(norm = norm)
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
        data_rgb = M2RGB(mxp,myp,mzp)
        ls = mcolors.LightSource()
        plt.imshow(np.transpose(data_rgb), origin="lower", interpolation="bicubic")
    end

    plt.xticks([])
    plt.yticks([])

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"m")
    if !isdir(path)
        mkpath(path)
    end
    plt.savefig(joinpath(path,basename(fname))*"_$style.png",dpi=300)
    plt.close()
end