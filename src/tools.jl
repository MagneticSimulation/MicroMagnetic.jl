using FFTW

function fft_m(spin, nx, ny, nz, dx, dy, dz, axis)
    if !(axis in ['z', 'x', 'y'])
        error("axis should be one of 'x', 'y' and 'z'!!!")
    end

    m = reshape(spin,(3, nx, ny, nz))
    c = Int(axis) - Int('x') + 1
    fft_m = fftshift(fft(m[c,:,:,:]))

    idx = floor(Int,nx/2+1)
    idy = floor(Int,ny/2+1)
    idz = floor(Int,nz/2+1)

    kx = fftshift(fftfreq(nx, d=dx)*2*pi)
    ky = fftshift(fftfreq(ny, d=dy)*2*pi)
    kz = fftshift(fftfreq(nz, d=dz)*2*pi)

    if c == 3
        Intensity = zeros(nx,ny)
        for i =1:nx, j = 1:ny
            Intensity[i,j] = sum((abs.(fft_m[i,j,:])).^2)/nz
        end
        Intensity[idx,idy] = 0
        return kx, ky, Intensity
    elseif c == 2
        Intensity = zeros(nz,nx)
        for k = 1:nz, i = 1:nx
            Intensity[k,i] = sum((abs.(fft_m[k,:,i])).^2)/ny
        end
        Intensity[idz,idx] = 0
        return kz, kx, Intensity
    elseif c == 1
        Intensity = zeros(ny,nz)
        for j=1:ny, k =1:nz
            Intensity[j,k] = sum((abs.(fft_m[:,j,k])).^2)/nx
        end
        Intensity[idy,idz] = 0
        return ky, kz, Intensity
    end

    return nothing
end

function fft_m(ovf_name; axis='z')

    ovf = read_ovf(ovf_name)
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize
    nxyz = nx*ny*nz
    spin = ovf.data

    return fft_m(spin, nx, ny, nz, dx, dy, dz, axis)
end
