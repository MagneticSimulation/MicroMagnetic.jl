using FFTW
using PyCall
"""
    OVF2XRAY(fname; xlim="none", ylim="none", axis=ez, zero_padding_size=-1)

Plot the simulated synchrotron radiation or neutron scattering contrast map.

And return the intensity array.

The plot boundary can be set by using "xlim = [-1e9,1e9]" "ylim = [-1e9,1e9]"
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
    nxyz = nx*ny*nz
    spin = ovf.data

    mxp, myp, mzp = sum_ovf(ovf, axis=ez)
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

## for test?
function fft_m_2d_direct(data, dx, dy, kxs, kys)

    nx, ny = size(data)
    fft_m = zeros(Complex{Float64}, length(kxs), length(kys))
    for I in 1:length(kxs), J in 1:length(kys)
        kx, ky = kxs[I], kys[J]
        for i = 1:nx, j=1:ny
            fft_m[I,J] += data[i,j]*exp(-1im*(kx*i*dx+ky*j*dy))
        end
    end
    return (abs.(fft_m)).^2
end

function fft_m(ovf_name, kxs, kys; axis='z')

    ovf = read_ovf(ovf_name)
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize
    nxyz = nx*ny*nz
    spin = ovf.data

    m = reshape(spin,(3, nx, ny, nz))
    c = Int(axis) - Int('x') + 1

    Intensity = zeros(length(kxs),length(kys))
    if c == 3
        for k = 1:nz
            I = fft_m_2d_direct(m[3,:,:,k], dx, dy, kxs, kys)
            Intensity .+= I
        end
        return Intensity./nz
    elseif c == 2
        println("will be added later!")
        return nothing
    elseif c == 1
        for i = 1:nx
            I = fft_m_2d_direct(m[1,i,:,:], dy, dz, kxs, kys)
            Intensity .+= I
        end
        return Intensity./nx
    end

    return fft_m(spin, nx, ny, nz, dx, dy, dz, axis)
end