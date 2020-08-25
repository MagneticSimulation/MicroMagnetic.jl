using FFTW
using PyCall

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

"""
        OVF2LTEM(fname;df=200, Ms=1e5, V=300, V0=-26, alpha=1e-5, Nx=-1, Ny=-1,beta=0,gamma=0)

Plot the simulated LTEM contrast from a .ovf file.

df: defocus length in um.

Ms: magnetization

Nx,Ny: output size

beta: Rotation angle with minus-y-axis (0,-1,0) in radian

gamma: Rotation angle with x-axis in radian

+:anti-clockwise -:clockwise

V: the Accelerating voltage in Kv

V0: the mean inner potential (MIP)
"""
function OVF2LTEM(fname;df=200, Ms=1e5, V=300, V0=-26, alpha=1e-5, Nx=-1, Ny=-1,beta=0,gamma=0)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")
    colorsys = pyimport("colorsys")
    ag = pyimport("mpl_toolkits.axes_grid1")

    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"LTEM")

    if !isdir(path)
        mkpath(path)
    end

    phase, intensity = LTEM(fname; V=V, Ms=1e5, V0=V0, df=df, alpha=alpha, Nx=Nx, Ny=Ny,beta=beta,gamma=gamma)
    fig,ax = plt.subplots(1,2)
    ax[1].imshow(transpose(phase),cmap=mpl.cm.gray,origin="lower") 
    ax[1].set_title("phase")
    ax[2].imshow(transpose(intensity),cmap=mpl.cm.gray,origin="lower")
    ax[2].set_title("intensity")
    plt.savefig(joinpath(path,basename(fname)*"_LTEM.png"),dpi=300)
    plt.close()

    return phase, intensity
end

#Ref: "MALTS: A tool to simulate Lorentz Transmission Electron Microscopy from micromagnetic simulations" by Stephanie K. Walton
#df in um
#V is the Accelerating voltage, in Kv
#V0 is the mean inner potential (MIP)
#alpha: beam divergence angle
#axis="x" is the normal direction
#beta is the angle between electron beam and the normal  axis, in radian
function LTEM(fname; V=300, Ms=1e5, V0=-26, df=1600, alpha=1e-5, Nx=-1, Ny=-1,beta=0,gamma=0)
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

    Nx = Nx > nx ? Nx : nx
    Ny = Ny > ny ? Ny : ny

    sum_mx,sum_my,sum_mz = Make_Projection(ovf,beta=beta,gamma=gamma,Nx=Nx,Ny=Ny)
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

    #if i,j within material boundary, add electric phase

    for i=1:Nx,j=1:Ny
        if abs(smx[i,j]) + abs(smy[i,j]) +abs(smz[i,j]) > 1e-5
            phi[i,j] += phi_E
        end
    end
    phi = mod.(phi,2*pi)

    fg = fft(exp.(1im.*phi))

    intensity = (abs.(ifft(fg.*E.*T))).^2

    #return phi_M, intensity
    return phi, intensity
end