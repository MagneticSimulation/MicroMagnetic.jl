using FFTW
#using PyCall

function fft_m_2d(data, dx, dy, zero_padding_size)

    N = zero_padding_size
    nx, ny = size(data)
    Nx = N > nx ? N : nx
    Ny = N > ny ? N : ny

    new_data = zeros(Nx, Ny)
    for i=1:nx, j=1:ny
        new_data[i,j] = data[i,j]
    end

    fft_m = fftshift(fft(new_data))

    idx = floor(Int,nx/2+1)
    idy = floor(Int,ny/2+1)

    kx = fftshift(fftfreq(Nx, d=dx)*2*pi)
    ky = fftshift(fftfreq(Ny, d=dy)*2*pi)

    Intensity = (abs.(fft_m)).^2

    return kx, ky, Intensity
end


function fft_m2(spin, nx, ny, nz, dx, dy, dz, axis; zero_padding_size=-1)
    if !(axis in [ex,ey,ez])
        error("axis should be one of 'ex', 'ey' and 'ez'!!!")
    end

    m = reshape(spin,(3, nx, ny, nz))

    idx = floor(Int,nx/2+1)
    idy = floor(Int,ny/2+1)
    idz = floor(Int,nz/2+1)

    N = zero_padding_size
    Nx = N > nx ? N : nx
    Ny = N > ny ? N : ny
    Nz = N > nz ? N : nz

    kx = fftshift(fftfreq(Nx, d=dx)*2*pi)
    ky = fftshift(fftfreq(Ny, d=dy)*2*pi)
    kz = fftshift(fftfreq(Nz, d=dz)*2*pi)


    if axis == ez
        Intensity = zeros(Nx, Ny)
        for k = 1:nz
            mx,my,mz = m[1,:,:,k],m[2,:,:,k],m[3,:,:,k]
            pmx,pmy,pmz = Normal_Projection_z(mx,my,mz)
            kx, ky, I = fft_m_2d(pmz, dx, dy, zero_padding_size)
            Intensity .+= I
        end
        return kx, ky, Intensity./nz
    elseif axis == ey
        Intensity = zeros(Nx, Nz)
        for j = 1:ny
            mx,my,mz = m[1,:,j,:],m[2,:,j,:],m[3,:,j,:]
            pmx,pmy,pmz = Normal_Projection_y(mx,my,mz)
            kx, kz, I = fft_m_2d(pmz, dx, dz, zero_padding_size)
            Intensity .+= I
        end
        return kx, kz, Intensity./ny
    elseif axis == ex
        Intensity = zeros(Ny, Nz)
        for i = 1:nx
            mx,my,mz = m[1,i,:,:],m[2,i,:,:],m[3,i,:,:]
            pmx,pmy,pmz = Normal_Projection_x(mx,my,mz)
            ky, kz, I = fft_m_2d(pmz, dy, dz, zero_padding_size)
            Intensity .+= I
        end
        #Intensity[idy,idz] = 0
        return ky, kz, Intensity./nx
    end

    return nothing
end
"""
    OVF2XRAY(fname; xlim="none", ylim="none", axis=ez, zero_padding_size=-1)

Plot the simulated synchrotron radiation or neutron scattering contrast map.

And return the intensity array.

The plot boundary can be set by using "xlim = [-1e9,1e9]" "ylim = [-1e9,1e9]"
"""
function OVF2XRAY(fname; xlim="none", ylim="none", axis=ez, zero_padding_size=-1)
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

    kx, ky, Intensity = fft_m2(spin, nx, ny, nz, dx, dy, dz, axis, zero_padding_size=zero_padding_size)

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
