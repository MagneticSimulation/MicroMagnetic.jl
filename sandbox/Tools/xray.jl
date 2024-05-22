using FFTW
using PyCall
using MicroMagnetic

"""
    OVF2XRAY(fname; axis=ez, N=-1)

Plot the simulated synchrotron radiation or neutron scattering contrast map.

And return the intensity array.

The maximum shown wavelength can be limited by using "xlim=1e-9"
"""
function OVF2XRAY(fname; xlim="none", ylim="none", axis=ez, N=-1)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    plt = pyimport("matplotlib.pyplot")

    ovf = read_ovf(fname)
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize
    n_total = nx*ny*nz
    spin = ovf.data

    mxp, myp, mzp = sum_ovf(ovf, axis=axis)
    lx, ly = size(mzp)
    if N > 0
        mzp = np.pad(mzp, ((0,N-lx),(0,N-ly)),"constant")
        mzp = np.roll(mzp, (-floor(int, nx/2), -floor(int, ny/2)))
    end
    Nx, Ny = size(mzp)
    if axis == ez
        Dx, Dy = dx, dy
    elseif axis == ex
        Dx, Dy = dy, dz
    elseif axis == ey
        Dx, Dy = dx, dz
    end
    Dx, Dy = Dx * Nx, Dy * Ny

    Intensity = np.real(np.fft.fft2(mzp))
    Intensity = fftshift(Intensity)
    kx = fftshift(fftfreq(Nx, d=Dx)*2*pi)
    ky = fftshift(fftfreq(Ny, d=Dy)*2*pi)

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"XRAY")

    if !isdir(path)
        mkpath(path)
    end

    plt.imshow(np.transpose(Intensity), interpolation = "gaussian", origin="lower",
                            extent=(kx[1], kx[end], ky[1], ky[end]), cmap=mpl.cm.PuBu_r)
    if xlim != "none"
        plt.xlim(xlim)
    end
    if ylim != "none"
        plt.ylim(xlim)
    end
    plt.xlabel("kx")
    plt.ylabel("ky")
    plt.grid()
    cbar = plt.colorbar()

    plt.savefig(joinpath(path,basename(fname)*"_XRAY.png"),dpi=300)
    plt.close()

    return Intensity
end
