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


function compute_electric_phase(V, V0, dz, nz, beta)
    C = 299792458.0
    E0 = m_e*C^2
    Ek = V*c_e
    P =  sqrt(Ek^2+2*Ek*E0)/C
    lambda = 2*pi*h_bar/P #lambda in m
    CE = (2*pi/(lambda*V))*(Ek+E0)/(Ek+2*E0)
    phi_E = CE*V0*(dz*nz)/cos(beta)
    return lambda, phi_E
end

function compute_magnetic_phase_direct(mx, my, dx, dy, Lz, Ms, i0, j0)  #slow, for testing purpose
    Nx, Ny = size(mx)
    mu0 = 4*pi*1e-7
    Phi0 = 2.067833e-15

    phi = 0
    for i=1:Nx, j=1:Ny
        x = (i0-i)*dx
        y = (j0-j)*dy
        r = sqrt(x^2+y^2)
        if  r > 0.001*min(dx, dy)
            phi += (y*mx[i,j] -  x*my[i,j])/r^2
        end
    end
    return -phi*mu0*Ms*Lz/(2*Phi0)*dx*dy
end


function compute_magnetic_phase_fft(mx, my, dx, dy, Lz, Ms; Nx=-1, Ny=-1)
    mu0 = 4*pi*1e-7
    Phi0 = 2.067833e-15
    nx, ny =  size(mx)

    Nx = Nx < 0 ? nx : Nx
    Ny = Ny < 0 ? ny : Ny

    Lx = Nx + nx - 1
    Ly = Ny + ny - 1
    new_mx = zeros(Lx, Ly)
    new_my = zeros(Lx, Ly)
    for i=1:nx, j=1:ny
        new_mx[i,j] = mx[i,j]
        new_my[i,j] = my[i,j]
    end

    F(x::Number, y::Number) = y == 0 ? 0 : y/(x^2+y^2)
    G(x::Number, y::Number) = x == 0 ? 0 : x/(x^2+y^2)

    fs = zeros(Lx, Ly)
    gs = zeros(Lx, Ly)
    for i=1:Lx, j=1:Ly
        fs[i,j] = F((-(Lx-1)/2+i-1)*dx, (-(Ly-1)/2+j-1)*dy)
        gs[i,j] = G((-(Lx-1)/2+i-1)*dx, (-(Ly-1)/2+j-1)*dy)
    end

    fft_mx = fft(new_mx)
    fft_my = fft(new_my)
    fft_fs = fft(fs)
    fft_gs = fft(gs)

    fft_phi = zeros(Complex{Float64}, (Lx,Ly))
    for i=1:Lx, j=1:Ly
        fft_phi[i,j] = fft_mx[i,j]*fft_fs[i,j] - fft_my[i,j]*fft_gs[i,j]
    end

    phi_all = real.(ifft(fft_phi))
    phi = zeros(Nx, Ny)
    for i=1:Nx, j=1:Ny
        phi[i,j] = phi_all[i+nx-1, j+ny-1]
    end

    return -phi*mu0*Ms*Lz/(2*Phi0)*dx*dy
end

#Ref: "MALTS: A tool to simulate Lorentz Transmission Electron Microscopy from micromagnetic simulations" by Stephanie K. Walton
#df in um
#V is the Accelerating voltage, in Kv
#V0 is the mean inner potential (MIP)
#alpha: beam divergence angle
#axis="x" is the normal direction
#beta is the angle between electron beam and the normal axis, in radian
function LTEM(ovf_name; V=300, Ms=1e5, V0=-26, df=1600, alpha=1e-5, zero_padding_size=512,axis="z",beta=0)
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

    mx1,my1,mz1 = m_average(ovf,axis=axis)
    mx = mx1
    my = cos(beta)*my1-sin(beta)*mz1
    if axis == "x"
        nx,dx=ny,dy
        ny,dy=nz,dz
    elseif axis == "y"
        ny,dy=nz,dz
    end

    dx=dx*cos(beta)

    lambda, phi_E = compute_electric_phase(1000*V, V0, dz, nz, beta)

    N = zero_padding_size
    Nx = N > nx ? N : nx
    Ny = N > ny ? N : ny

    #put data on the center of zero padding square
    new_mx = zeros(Nx, Ny)
    new_my = zeros(Nx, Ny)
    for i=1:nx, j=1:ny
        new_mx[Int(floor(Nx/2-nx/2+i)),Int(floor(Ny/2-ny/2+j))] = mx[i,j]
        new_my[Int(floor(Nx/2-nx/2+i)),Int(floor(Ny/2-ny/2+j))] = my[i,j]
    end

    fft_mx = fft(new_mx)
    fft_my = fft(new_my)

    kx = fftfreq(Nx, d=dx)
    ky = fftfreq(Ny, d=dy)

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
        T[i,j] = exp(-pi*1im*(df*lambda*k2))
        E[i,j] = exp(-(pi*alpha*df*k)^2)
        if k2 > 0
            Phi_M[i,j] = 1im*(c_e/h_bar)*mu_0*Ms*(nz*dz)*(fft_mx_ky[i,j]-fft_my_kx[i,j])/k2
        end
    end
    Phi_M[1,1]=0

    phi_M = 1/cos(beta)*real.(ifft(Phi_M))

    phi = phi_M

    for i=1:nx,j=1:ny
        phi[Int(floor(Nx/2-nx/2+i)),Int(floor(Ny/2-ny/2+j))] += phi_E
    end
    phi = mod.(phi,2*pi)

    fg = fft(exp.(1im.*phi))

    intensity = (abs.(ifft(fg.*E.*T))).^2;

    local_phi = zeros(nx,ny)
    local_intensity = zeros(nx,ny)

    for i =1:nx,j=1:ny
        local_phi[i,j] = phi_M[Int(floor(Nx/2-nx/2+i)),Int(floor(Ny/2-ny/2+j))]
        local_intensity[i,j] = intensity[Int(floor(Nx/2-nx/2+i)),Int(floor(Ny/2-ny/2+j))]
    end

    return local_phi, local_intensity
end

function m_average(m::Array{T,1},nx::Int,ny::Int,nz::Int;axis::String="z") where T<:AbstractFloat ##axis can only chosen from "x" "y" "z"

    if length(m) != 3*nx*ny*nz
        println("Length doesn't match!")
        return nothing
    end

    b = reshape(m,(3,nx,ny,nz))
    N = 0
    if axis =="x"
        mx,my,mz = zeros(ny,nz),zeros(ny,nz),zeros(ny,nz)
        for i = 1:nx
            mx .+= b[2,i,:,:]
            my .+= b[3,i,:,:]
            mz .+= b[1,i,:,:]
        end
        N = nx
    elseif axis == "y"
        mx,my,mz = zeros(nx,nz),zeros(nx,nz),zeros(nx,nz)
        for j = 1:ny
            mx .-= b[1,:,j,:]
            my .+= b[3,:,j,:]
            mz .+= b[2,:,j,:]
        end
        N = ny
    elseif axis == "z"
        mx,my,mz = zeros(nx,ny),zeros(nx,ny),zeros(nx,ny)
        for k = 1:nz
            mx .+= b[1,:,:,k]
            my .+= b[2,:,:,k]
            mz .+= b[3,:,:,k]
        end
        N = nz
    end

    return mx/N,my/N,mz/N
end

function m_average(ovf::OVF2;axis::String="z")
    m = ovf.data
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes

    return m_average(m,nx,ny,nz,axis=axis)
end
