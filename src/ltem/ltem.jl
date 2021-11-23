using NPZ
using FFTW
using Printf

function vector_potential_constant(Ms::Float64, d::Float64)
    mu0 = 4*pi*1e-7
    Phi0 = 2.067833e-15
    return mu0*Ms/(2*Phi0*d)
end

function cross_product(v1::Array{T,4},v2::Array{T,4}) where {T<:Number}
    (three, nx, ny, nz) = size(v1)
    v3 = zeros(T, size(v1))

    for i = 1:nx, j=1:ny, k=1:nz
        v3[1,i,j,k] = cross_x(v1[1,i,j,k], v1[2,i,j,k], v1[3,i,j,k], v2[1,i,j,k], v2[2,i,j,k], v2[3,i,j,k])
        v3[2,i,j,k] = cross_y(v1[1,i,j,k], v1[2,i,j,k], v1[3,i,j,k], v2[1,i,j,k], v2[2,i,j,k], v2[3,i,j,k])
        v3[3,i,j,k] = cross_z(v1[1,i,j,k], v1[2,i,j,k], v1[3,i,j,k], v2[1,i,j,k], v2[2,i,j,k], v2[3,i,j,k])
    end
    return v3
end

"""
    get_magnetic_phase(m::Array{T,4}; alphas::Union{Array, Tuple}=[0],betas::Union{Array, Tuple}=[],N=256) where {T<:Number}

Get the LTEM phase shift from a magnetization array.

Input: (3, nx, ny, nz)

alphas: Tilt angles with axis X

betas: Tilt angles with axis Y

N: zeros padding size

"""

function get_magnetic_phase(m::Array{T,4}; alphas::Union{Array, Tuple}=[0],betas::Union{Array, Tuple}=[],N=256) where {T<:Number}

    a = compute_vector_potential(m, N=N)

    phase = vector_field_projection(a, alphas, "x")
    for i = 1:length(alphas)
        npzwrite(@sprintf("X_%g.npy",alphas[i]),phase[:,:,i])
    end

    phase = vector_field_projection(a, betas, "x")
    for i = 1:length(betas)
        npzwrite(@sprintf("X_%g.npy",betas[i]),phase[:,:,i])
    end
end

"""
    compute_vector_potential(m::Array{T,4}; N=256)  where {T<:Number}

Compute the magnetic vector potential from a magnetization array.

Input: (3, nx, ny, nz)

output: (3, N, N, N)

N: zeros padding size

"""

function compute_vector_potential(m::Array{T,4}; N=256)  where {T<:Number}
    m = vector_padding(m,N,N,N)
    m = Array{Float32}(m)
    m = ifftshift(m, (2,3,4))
    mk = fft(m, (2,3,4))
    k = get_kernel_k(N)
    vector_potential_k = cross_product(mk,k)
    vector_potential = ifft(vector_potential_k, (2,3,4))
    return fftshift(real.(vector_potential), (2,3,4))
end

function get_kernel_k(N::Int)
    r = FFTW.fftfreq(N)
    kernel_k = zeros(ComplexF32, (3,N,N,N))
    for i = 1:N
        x2 = r[i]^2
        for j = 1:N
            y2 = r[j]^2
            for k = 1:N
                z2 = r[k]^2
                r2 = x2+y2+z2

                kernel_k[1,i,j,k] = -1im*r[i]/r2
                kernel_k[2,i,j,k] = -1im*r[j]/r2
                kernel_k[3,i,j,k] = -1im*r[k]/r2
            end
        end
    end
    kernel_k[1,1,1,1] = 0
    kernel_k[2,1,1,1] = 0
    kernel_k[3,1,1,1] = 0
    return kernel_k        
end

function compute_electric_phase(V, V0, Lz, beta)
    C = 299792458.0
    E0 = m_e*C^2
    Ek = V*c_e
    P =  sqrt(Ek^2+2*Ek*E0)/C
    lambda = 2*pi*h_bar/P #lambda in m
    CE = (2*pi/(lambda*V))*(Ek+E0)/(Ek+2*E0)
    phi_E = CE*V0*Lz/cos(deg2rad(beta))
    return lambda, phi_E
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
#=function OVF2LTEM(fname;df=200, Ms=1e5, V=300, V0=-26, alpha=1e-5, N=-1,tilt_axis="y",tilt_angle=0)
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

    phase, intensity = LTEM(fname; V=V, Ms=Ms, V0=V0, df=df, alpha=alpha, N=N, tilt_axis=tilt_axis, tilt_angle=tilt_angle)
    fig,ax = plt.subplots(1,2)
    ax[1].imshow(np.transpose(fftshift(phase)),cmap=mpl.cm.gray,origin="lower")
    ax[1].set_title("phase")
    ax[2].imshow(np.transpose(fftshift(intensity)),cmap=mpl.cm.gray,origin="lower")
    ax[2].set_title("intensity")
    plt.savefig(joinpath(path,basename(fname)*"_LTEM.png"),dpi=300)
    plt.close()

    return phase, intensity
end=#

#Ref: "MALTS: A tool to simulate Lorentz Transmission Electron Microscopy from micromagnetic simulations" by Stephanie K. Walton
#df in um
#V is the Accelerating voltage, in Kv
#V0 is the mean inner potential (MIP)
#alpha: beam divergence angle
#axis="x" is the normal direction
#beta is the angle between electron beam and the normal  axis, in radian
#=function LTEM(fname; V=300, Ms=1e5, V0=-26, df=1600, alpha=1e-5, N=128,tilt_axis="y",tilt_angle=0)
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

    if nx > N || ny > N
        @error("Please add zero padding!")
    end

    mp = radon_transform_ovf(ovf, tilt_angle, tilt_axis,N=N)
    sum_mx, sum_my, sum_mz = mp[1,:,:], mp[2,:,:], mp[3,:,:]
    smx, smy, smz = sum_mx/nz, sum_my/nz, sum_mz/nz

    phi_M = compute_magnetic_phase_fft(smx, smy, dx, dy, dz, Ms)

    lambda, phi_E = compute_electric_phase(1000*V, V0, dz*nz, tilt_angle)

    #put data on the center of zero padding square


    kx = fftfreq(N, d=dx)
    ky = fftfreq(N, d=dy)


    Phi_M = zeros(Complex{Float64}, (N,N))
    T = zeros(Complex{Float64}, (N,N))
    E = zeros(N,N)
    df = df*1e-6

    for i=1:N, j=1:N
        k2 = kx[i]^2 + ky[j]^2
        k = sqrt(k2)
        T[i,j] = exp(-pi*1im*(df*lambda*k2))
        E[i,j] = exp(-(pi*alpha*df*k)^2)
    end


    #if i,j within material boundary, add electric phase

    #=for i=1:N,j=1:N
        if abs(smx[i,j]) + abs(smy[i,j]) +abs(smz[i,j]) > 1e-5
            phi[i,j] += phi_E
        end
    end=#
    #phi = mod.(phi,2*pi)

    fg = fft(exp.(1im.*phi))

    intensity = (abs.(ifft(fg.*E.*T))).^2

    #return phi_M, intensity
    return phi, intensity
end=#
