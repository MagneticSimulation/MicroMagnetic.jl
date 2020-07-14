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
        Intensity = zeros(Nx, Nz)
        for j = 1:ny
            kx, kz, I = fft_m_2d(m[2,:,j,:], dx, dz, zero_padding_size)
            Intensity .+= I
        end
        return kx, kz, Intensity./ny
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

function compute_magnetic_phase_direct(mx, my, dx, dy, dz, Ms, i0, j0)  #slow, for testing purpose
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
    return -phi*mu0*Ms*dz/(2*Phi0)*dx*dy
end


function compute_magnetic_phase_fft(mx, my, dx, dy, dz, Ms; Nx=-1, Ny=-1)
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

    return -phi*mu0*Ms*dz/(2*Phi0)*dx*dy
end

#Ref: "MALTS: A tool to simulate Lorentz Transmission Electron Microscopy from micromagnetic simulations" by Stephanie K. Walton
#df in um
#V is the Accelerating voltage, in Kv
#V0 is the mean inner potential (MIP)
#alpha: beam divergence angle
#axis="x" is the normal direction
#beta is the angle between electron beam and the normal  axis, in radian
function LTEM(fname; V=300, Ms=1e5, V0=-26, df=1600, alpha=1e-5, Nx=512, Ny=512,beta=0,gamma=0,ItpNum=0)
    ovf = read_ovf(fname)
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize
    nxyz = nx*ny*nz
    spin = ovf.data
    m = reshape(spin,(3, nx, ny, nz))

    sum_mx,sum_my,sum_mz = Make_Projection(ovf,beta=beta,gamma=gamma,Nx=Nx,Ny=Ny,ItpNum=ItpNum)
    smx,smy,smz = sum_mx/nz,sum_my/nz,sum_mz/nz

    #dx=dx*cos(beta)

    lambda, phi_E = compute_electric_phase(1000*V, V0, dz, nz, beta)

    #put data on the center of zero padding square

    fft_mx = fft(smx)
    fft_my = fft(smy)

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

    phi_M = real.(ifft(Phi_M))

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
    #return local_phi, intensity
end

"""Make_Projection: get the projection magnetization.
nx ny nz is the shape of m, Nx Ny is the projection size

First step: rotate beta(rad) with axis (0,-1,0)
Second step : rotate gamma(rad) with axis (1,0,0)
positive degree for anti-clockwise,negative for clockwise

ItpNum:Interpolation within an unit length
for example:ItpNum=2,insert 3*3*3-8 unit cells from [i,j,k] to [i+1,j+1,k+1] for every cell
Using ItpNum can improve the accuracy of projection operation at a cost of calculation speed
(not suggested in ordinary)

```julia
Make_Projection(m,nx,ny,nz,Nx=128,Ny=128,beta=pi/3)
```
"""
function Make_Projection(m::Array{T,1},nx,ny,nz; Nx::Int=-1,Ny::Int=-1,
                        beta::Number=0,gamma::Number=0,ItpNum::Int=0) where T<:AbstractFloat

    Nx = Nx > 0 ? Nx : nx
    Ny = Ny > 0 ? Ny : ny

    b = reshape(m,(3,nx,ny,nz))
    mx,my,mz = b[1,:,:,:],b[2,:,:,:],b[3,:,:,:]

    mxp,myp,mzp = Make_Projection(mx,my,mz,Nx=Nx,Ny=Ny,beta=beta,gamma=gamma,ItpNum=ItpNum)
    return mxp,myp,mzp
end

function trilinear_interpolation(M,x,y,z)
    (nx,ny,nz)=size(M)
    xf,yf,zf = Int(floor(x)),Int(floor(y)),Int(floor(z))
    if x==nx || y==ny || z==nz
        return M[xf,yf,zf]
    end 
    xd,yd,zd = x-xf,y-yf,z-zf
    xp,yp,zp = 1-xd,1-yd,1-zd

    v = M[xf,yf,zf]*xp*yp*zp + M[xf+1,yf+1,zf+1]*xd*yd*zd +
        M[xf,yf,zf+1]*xp*yp*zd + M[xf,yf+1,zf]*xp*yd*zp + M[xf+1,yf,zf]*xd*yp*zp+
        M[xf+1,yf+1,zf]*xp*yd*zp + M[xf+1,yf,zf+1]*xd*yp*zd + M[xf,yf+1,zf+1]*xp*yd*zd

    return v
end

function Make_Projection(mx::Array{Float64,3},my::Array{Float64,3},mz::Array{Float64,3};
                          Nx::Int=128,Ny::Int=128,beta::Number=0,gamma::Number=0,ItpNum::Int=0)

    (nx,ny,nz) = size(mx)
    mxp,myp,mzp = zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)
    eps0 = 1/(ItpNum+1)
    eps3 = eps0*eps0*eps0
    for i=1:nx, j= 1:ny, k=1:nz
        for eps_i = 0:ItpNum, eps_j = 0:ItpNum, eps_k = 0:ItpNum
            x0,y0,z0 = i+eps_i*eps0, j+eps_j*eps0, k+eps_k*eps0
            if x0>nx+1e-6 || y0>ny+1e-6 || z0>nz+1e-6
                continue
            end

            if eps_i==0 && eps_j==0 && eps_k==0
                itp_mx,itp_my,itp_mz = mx[i,j,k],my[i,j,k],mz[i,j,k]
            else
                itp_mx = trilinear_interpolation(mx,x0,y0,z0)
                itp_my = trilinear_interpolation(my,x0,y0,z0)
                itp_mz = trilinear_interpolation(mz,x0,y0,z0)
            end

            x1,y1,z1= x0-(nx+1)/2,y0-(ny+1)/2,z0-(nz+1)/2
            axis1 = [0,-1.0,0]
            x2,y2,z2 = rotation_operator([x1,y1,z1],axis1,float(beta))

            axis2 = rotation_operator([1.0,0,0],axis1,float(beta))
            x3,y3,z3 = rotation_operator([x2,y2,z2],axis2,float(gamma))
            x,y = x3+(Nx+1)/2,y3+(Ny+1)/2
            linear_separate(mxp,myp,mzp,itp_mx,itp_my,itp_mz,x,y,Nx,Ny,beta,gamma)
        end
    end
    return mxp*eps3,myp*eps3,mzp*eps3
end

function linear_separate(mxp,myp,mzp,mx,my,mz,x,y,Nx,Ny,beta,gamma)
    cb,sb,cg,sg = cos(beta),sin(beta),cos(gamma),sin(gamma)
    if 1 <= x <= Nx && 1 <= y <= Ny
        xf,yf = Int(floor(x)),Int(floor(y))
        xd,yd = x-xf,y-yf
        xp,yp = 1-xd,1-yd
        mx1,my1,mz1 = mx*cb-mz*sb, my, mz*cb+mx*sb
        mxx,myy,mzz = mx1, my1*cg-mz1*sg, my1*sg+mz1*cg

        mxp[xf,yf] += mxx*xp*yp
        myp[xf,yf] += myy*xp*yp
        mzp[xf,yf] += mzz*xp*yp

        if yf+1 <= Ny
            mxp[xf,yf+1] += mxx*xp*yd
            myp[xf,yf+1] += myy*xp*yd
            mzp[xf,yf+1] += mzz*xp*yd
        end

        if xf+1 <= Nx
            mxp[xf+1,yf] += mxx*xd*yp
            myp[xf+1,yf] += myy*xd*yp
            mzp[xf+1,yf] += mzz*xd*yp
        end
        
        if yf+1 <= Ny && xf+1 <= Nx
            mxp[xf+1,yf+1] += mxx*xd*yd
            myp[xf+1,yf+1] += myy*xd*yd
            mzp[xf+1,yf+1] += mzz*xd*yd
        end
    end
end

function Make_Projection(ovf::OVF2;beta=0,gamma=0,Nx=256,Ny=256,ItpNum=0)
    m = ovf.data
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes

    return Make_Projection(m,nx,ny,nz,beta=beta,gamma=gamma,Nx=Nx,Ny=Ny,ItpNum=ItpNum)
end


