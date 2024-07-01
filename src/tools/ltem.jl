using FFTW
using Printf

export compute_magnetic_phase

function compute_magnetic_phase_direct(mx, my, dx, dy, Lz, Ms, i0, j0)  #slow, for testing purpose
    Nx, Ny = size(mx)
    mu0 = 4 * pi * 1e-7
    Phi0 = 2.067833e-15

    phi = 0
    for i in 1:Nx, j in 1:Ny
        x = (i0 - i) * dx
        y = (j0 - j) * dy
        r = sqrt(x^2 + y^2)
        if r > 0.001 * min(dx, dy)
            phi += (y * mx[i, j] - x * my[i, j]) / r^2
        end
    end
    return -phi * mu0 * Ms * Lz / (2 * Phi0) * dx * dy
end

function compute_magnetic_phase_fft(mx, my, dx, dy, Lz, Ms; Nx=-1, Ny=-1)
    mu0 = 4 * pi * 1e-7
    Phi0 = 2.067833e-15
    nx, ny = size(mx)

    Nx = Nx < 0 ? nx : Nx
    Ny = Ny < 0 ? ny : Ny

    Lx = Nx + nx - 1
    Ly = Ny + ny - 1
    new_mx = zeros(Lx, Ly)
    new_my = zeros(Lx, Ly)
    for i in 1:nx, j in 1:ny
        new_mx[i, j] = mx[i, j]
        new_my[i, j] = my[i, j]
    end

    F(x::Number, y::Number) = y == 0 ? 0 : y / (x^2 + y^2)
    G(x::Number, y::Number) = x == 0 ? 0 : x / (x^2 + y^2)

    fs = zeros(Lx, Ly)
    gs = zeros(Lx, Ly)
    for i in 1:Lx, j in 1:Ly
        fs[i, j] = F((-(Lx - 1) / 2 + i - 1) * dx, (-(Ly - 1) / 2 + j - 1) * dy)
        gs[i, j] = G((-(Lx - 1) / 2 + i - 1) * dx, (-(Ly - 1) / 2 + j - 1) * dy)
    end

    fft_mx = fft(new_mx)
    fft_my = fft(new_my)
    fft_fs = fft(fs)
    fft_gs = fft(gs)

    fft_phi = zeros(Complex{Float64}, (Lx, Ly))
    for i in 1:Lx, j in 1:Ly
        fft_phi[i, j] = fft_mx[i, j] * fft_fs[i, j] - fft_my[i, j] * fft_gs[i, j]
    end

    phi_all = real.(ifft(fft_phi))
    phi = zeros(Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        phi[i, j] = phi_all[i + nx - 1, j + ny - 1]
    end

    return -phi * mu0 * Ms * Lz / (2 * Phi0) * dx * dy
end

function average_m(m::AbstractArray{T,1}, nx::Int, ny::Int, nz::Int;
                   axis::Symbol=:z) where {T<:AbstractFloat}
    # Check if the length of m matches the expected length
    if length(m) != 3 * nx * ny * nz
        println("Length doesn't match!")
        return nothing
    end

    b = reshape(m, (3, nx, ny, nz))
    axis_to_dim = Dict(:x => 2, :y => 3, :z => 4)
    dim = get(axis_to_dim, axis, nothing)
    if dim === nothing
        println("Invalid axis!")
        return nothing
    end

    # Calculate the average for each component
    average_b = sum(b; dims=dim) / size(b, dim)

    # Reshape the result to a 2D matrix
    if axis == :x
        result = reshape(average_b, (3, ny, nz))
    elseif axis == :y
        result = reshape(average_b, (3, nx, nz))
    elseif axis == :z
        result = reshape(average_b, (3, nx, ny))
    end

    return result[1, :, :], result[2, :, :], result[3, :, :]
end

function compute_magnetic_phase(ovf::OVF2; Ms=1, axis::Symbol=:z, Nx=-1, Ny=-1)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    mx, my, mz = average_m(ovf.data, nx, ny, nz; axis=axis)
    dx, dy, dz = ovf.xstepsize, ovf.ystepsize, ovf.zstepsize
    Lz = dz * ovf.znodes

    return compute_magnetic_phase_fft(mx, my, dx, dy, Lz, Ms; Nx=Nx, Ny=Ny)
end

# FIXME: we need to check the correctness of this function
#Ref: "MALTS: A tool to simulate Lorentz Transmission Electron Microscopy from micromagnetic simulations" by Stephanie K. Walton
#df in um
#V is the Accelerating voltage, in Kv
#V0 is the mean inner potential (MIP)
#alpha: beam divergence angle
#axis="x" is the normal direction
#beta is the angle between electron beam and the normal axis, in radian
function LTEM(ovf_name; V=300, Ms=1e5, V0=-26, df=1600, alpha=1e-5, zero_padding_size=512,
              axis="z", beta=0)
    ovf = read_ovf(ovf_name)
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize
    n_total = nx * ny * nz
    spin = ovf.data
    m = reshape(spin, (3, nx, ny, nz))

    mx1, my1, mz1 = m_average(ovf; axis=axis)
    mx = mx1
    my = cos(beta) * my1 - sin(beta) * mz1
    if axis == "x"
        nx, dx = ny, dy
        ny, dy = nz, dz
    elseif axis == "y"
        ny, dy = nz, dz
    end

    dx = dx * cos(beta)

    lambda, phi_E = compute_electric_phase(1000 * V, V0, dz * nz, beta)

    N = zero_padding_size
    Nx = N > nx ? N : nx
    Ny = N > ny ? N : ny

    #put data on the center of zero padding square
    new_mx = zeros(Nx, Ny)
    new_my = zeros(Nx, Ny)
    for i in 1:nx, j in 1:ny
        new_mx[Int(floor(Nx / 2 - nx / 2 + i)), Int(floor(Ny / 2 - ny / 2 + j))] = mx[i, j]
        new_my[Int(floor(Nx / 2 - nx / 2 + i)), Int(floor(Ny / 2 - ny / 2 + j))] = my[i, j]
    end

    fft_mx = fft(new_mx)
    fft_my = fft(new_my)

    kx = fftfreq(Nx; d=dx)
    ky = fftfreq(Ny; d=dy)

    fft_mx_ky = zeros(Complex{Float64}, (Nx, Ny))
    fft_my_kx = zeros(Complex{Float64}, (Nx, Ny))
    for i in 1:Nx, j in 1:Ny
        fft_mx_ky[i, j] = fft_mx[i, j] * ky[j]
        fft_my_kx[i, j] = fft_my[i, j] * kx[i]
    end

    Phi_M = zeros(Complex{Float64}, (Nx, Ny))
    T = zeros(Complex{Float64}, (Nx, Ny))
    E = zeros(Nx, Ny)
    df = df * 1e-6

    for i in 1:Nx, j in 1:Ny
        k2 = kx[i]^2 + ky[j]^2
        k = sqrt(k2)
        T[i, j] = exp(-pi * 1im * (df * lambda * k2))
        E[i, j] = exp(-(pi * alpha * df * k)^2)
        if k2 > 0
            Phi_M[i, j] = 1im *
                          (c_e / h_bar) *
                          mu_0 *
                          Ms *
                          (nz * dz) *
                          (fft_mx_ky[i, j] - fft_my_kx[i, j]) / k2
        end
    end
    Phi_M[1, 1] = 0

    phi_M = 1 / cos(beta) * real.(ifft(Phi_M))

    phi = phi_M

    for i in 1:nx, j in 1:ny
        phi[Int(floor(Nx / 2 - nx / 2 + i)), Int(floor(Ny / 2 - ny / 2 + j))] += phi_E
    end
    phi = mod.(phi, 2 * pi)

    fg = fft(exp.(1im .* phi))

    intensity = (abs.(ifft(fg .* E .* T))) .^ 2

    local_phi = zeros(nx, ny)
    local_intensity = zeros(nx, ny)

    for i in 1:nx, j in 1:ny
        local_phi[i, j] = phi_M[Int(floor(Nx / 2 - nx / 2 + i)),
                                Int(floor(Ny / 2 - ny / 2 + j))]
        local_intensity[i, j] = intensity[Int(floor(Nx / 2 - nx / 2 + i)),
                                          Int(floor(Ny / 2 - ny / 2 + j))]
    end

    return local_phi, local_intensity
end

# FIXME: we need to check the correctness of this function
"""
    compute_magnetic_phase(m::Array{T,4}, theta::Real, axis::String;
N1::Int=-1, N2::Int=-1, Ms=1/mu_0, d::Real=1.0) where {T<:AbstractFloat}

Compute the magnetic phase from magnetization array.

Parameters
------------------------------
m: 4D array sized (3, nx, ny, nz). Magnetization array.
theta: Euler angle.
axis: rotation axis chosen from "X" or "Y"
N1: padding size when calculating projection
N2: padding size of Fourier transfrom kernel, by default 2*N1

Outputs
----------------------------
phi: 2D array sized (N,N). Magnetic phase.

"""
function compute_magnetic_phase(m::Array{T,4}, theta::Real, axis::String; N1::Int=-1,
                                N2::Int=-1, Ms=1 / mu_0,
                                d::Real=1.0) where {T<:AbstractFloat}
    (_, nx, ny, nz) = size(m)
    max_n = max(nx, ny, nz)
    N1 = N1 > max_n ? N1 : max_n
    N2 = N2 > N1 ? N2 : N1

    m = pad(m, (3, N1, N1, N1))
    mx, my, mz = radon3d_xyz(m, theta, axis)
    mx_pad, my_pad = pad(mx, (N2, N2)), pad(my, (N2, N2))

    mx_pad = ifftshift(mx_pad, (1, 2))
    mx_k = fft(mx_pad, (1, 2))

    my_pad = ifftshift(my_pad, (1, 2))
    my_k = fft(my_pad, (1, 2))

    ks = fftfreq(N2)
    phi_k = zeros(ComplexF32, N2, N2)
    for i in 1:N2, j in 1:N2
        k2 = ks[i]^2 + ks[j]^2
        # Akz need a -1, but projection also need a -1
        phi_k[i, j] = 1im * (mx_k[i, j] * ks[j] - my_k[i, j] * ks[i]) / k2
    end
    phi_k[1, 1] = 0.0
    Phi0 = 2.067833e-18
    coeff = mu_0 * Ms * d^2 * pi / Phi0
    phi = coeff .* real.(ifft(phi_k, (1, 2)))
    return fftshift(phi)
end

function compute_electric_phase(V, V0, Lz, beta)
    C = 299792458.0
    E0 = m_e * C^2
    Ek = V * c_e
    P = sqrt(Ek^2 + 2 * Ek * E0) / C
    lambda = 2 * pi * h_bar / P #lambda in m
    CE = (2 * pi / (lambda * V)) * (Ek + E0) / (Ek + 2 * E0)
    phi_E = CE * V0 * Lz / cos(deg2rad(beta))
    return lambda, phi_E
end
