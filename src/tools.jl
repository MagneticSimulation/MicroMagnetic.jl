using FFTW

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


function fft_m(spin, nx, ny, nz, dx, dy, dz, axis, zero_padding_size=-1)
    if !(axis in ['z', 'x', 'y'])
        error("axis should be one of 'x', 'y' and 'z'!!!")
    end

    m = reshape(spin,(3, nx, ny, nz))
    c = Int(axis) - Int('x') + 1

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


    if c == 3
        Intensity = zeros(Nx, Ny)
        for k = 1:nz
            kx, ky, I = fft_m_2d(m[3,:,:,k], dx, dy, zero_padding_size)
            Intensity .+= I
        end
        return kx, ky, Intensity./nz
    elseif c == 2
        println("will be added later!")
        return nothing
    elseif c == 1
        Intensity = zeros(Ny, Nz)
        for i = 1:nx
            ky, kz, I = fft_m_2d(m[1,i,:,:], dy, dz, zero_padding_size)
            Intensity .+= I
        end
        #Intensity[idy,idz] = 0
        return ky, kz, Intensity./nx
    end

    return nothing
end

function fft_m(ovf_name; axis='z', zero_padding_size=-1)

    ovf = read_ovf(ovf_name)
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize
    nxyz = nx*ny*nz
    spin = ovf.data

    return fft_m(spin, nx, ny, nz, dx, dy, dz, axis, zero_padding_size)
end


function compute_electric_phase(V, V0, dz, nz)
    C = 299792458.0
    E0 = m_e*C^2
    Ek = V*c_e
    P =  sqrt(Ek^2+2*Ek*E0)/C
    lambda = 2*pi*h_bar/P #lambda in m
    CE = (2*pi/(lambda*V))*(Ek+E0)/(Ek+2*E0)
    phi_E = CE*V0*(dz*nz)
    return lambda, phi_E
end


#df in um
#V is the Accelerating voltage, in Kv
#V0 is the mean inner potential (MIP)
#alpha: beam divergence angle
function LTEM(ovf_name; V=300, Ms=1e5, V0=-26, df=1600, Cs=0, alpha=1e-5, zero_padding_size=-1)
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

    lambda, phi_E = compute_electric_phase(1000*V, V0, dz, nz)

    mx = m[1, :, :, 1]
    my = m[2, :, :, 1]

    N = zero_padding_size
    Nx = N > nx ? N : nx
    Ny = N > ny ? N : ny

    new_mx = zeros(Nx, Ny)
    new_my = zeros(Nx, Ny)
    for i=1:nx, j=1:ny
        new_mx[i,j] = mx[i,j]
        new_my[i,j] = my[i,j]
    end

    fft_mx = fftshift(fft(new_mx))
    fft_my = fftshift(fft(new_my))

    kx = fftshift(fftfreq(Nx, d=dx)*2*pi)
    ky = fftshift(fftfreq(Ny, d=dy)*2*pi)

    fft_mx_ky = zeros(Complex{Float64}, (Nx,Ny))
    fft_my_kx = zeros(Complex{Float64}, (Nx,Ny))
    for i=1:Nx, j=1:Ny
        fft_mx_ky[i,j] = fft_mx[i,j]*ky[j]
        fft_my_kx[i,j] = fft_my[i,j]*kx[i]
    end

    Phi_M = zeros(Complex{Float64}, (Nx,Ny))
    T = zeros(Complex{Float64}, (Nx,Ny))
    E = zeros(Nx,Ny)
    df = df*1e-6
    for i=1:Nx, j=1:Ny
        k2 = kx[i]^2 + ky[j]^2
        k = sqrt(k2)
        T[i,j] = exp(-pi*1im*(0.5*Cs*lambda^3*k2-df*lambda*k2))
        E[i,j] = exp(-(pi^2)*(alpha^2)*(df*k + Cs*(lambda^2)*k^3)^2)
        if k2 > 0
            Phi_M[i,j] = 1im*(c_e/h_bar)*pi*mu_0*Ms*(nz*dz)*(fft_mx_ky[i,j]-fft_my_kx[i,j])/k2
        end
    end
    println(maximum(E))


    phi_M = real(ifft(Phi_M))
    println(maximum(phi_M), phi_E)
    phi = mod.(phi_M .+ phi_E, 2*pi)

    fg = exp.(1im.*phi)
    intensity = (abs.(ifft(fg.*T.*E))).^2;
    return phi_M, intensity
end
