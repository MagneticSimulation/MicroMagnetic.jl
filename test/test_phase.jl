using NuMag
using Printf
using Test

function F0(x::Float64, y::Float64)
    eps = 1e-10
    if abs(x)<eps && abs(y)<eps
        return 0
    end
    return x*log(x^2+y^2) - 2*x + 2*y*atan(x/y)
end

function phim_uniformly_magnetized_slab(x, y, mx, my, Lx, Ly, Lz, Ms)
    mu0 = 4*pi*1e-7
    Phi0 = 2.067833e-15
    a = F0(x-Lx/2, y-Ly/2) - F0(x+Lx/2, y-Ly/2) - F0(x-Lx/2, y+Ly/2) + F0(x+Lx/2, y+Ly/2)
    b = F0(y-Ly/2, x-Lx/2) - F0(y+Ly/2, x-Lx/2) - F0(y-Ly/2, x+Lx/2) + F0(y+Ly/2, x+Lx/2)
    return mu0*Ms*Lz/(4*Phi0)*(my*b-mx*a)*1e-18
end

function phase_in_theory()
    mx = cos(5/3*pi)
    my = sin(5/3*pi)
    Lx = 32
    Ly = 64
    Lz = 16
    Nx, Ny = 160, 160
    Ms = 1e5
    phi = zeros(Nx, Nx)
    for i = 1:Nx, j=1:Ny
        x = (i - div(Nx, 2)) - 0.5
        y = (j - div(Ny, 2)) - 0.5
        phi[i, j] = phim_uniformly_magnetized_slab(x, y, mx, my, Lx, Ly, Lz, Ms)
    end
    return phi
end

function phase_fft()
    mx = cos(5/3*pi)
    my = sin(5/3*pi)
    Lx = 32
    Ly = 64
    Lz = 16e-9
    Nx, Ny = 160, 160
    dx, dy = 1e-9, 1e-9

    mx_m = zeros(Lx, Ly)
    my_m = zeros(Lx, Ly)
    mx_m .= mx
    my_m .= my

    Ms = 1e5
    phi = NuMag.compute_magnetic_phase_fft(mx_m, my_m, dx, dy, Lz, Ms, Nx=Nx, Ny=Ny)
    return phi
end

function phase_direct()
    mx = cos(5/3*pi)
    my = sin(5/3*pi)
    Lx = 32
    Ly = 64
    Lz = 16e-9
    Nx, Ny = 160, 160
    dx, dy = 1e-9, 1e-9

    mx_m = zeros(Lx, Ly)
    my_m = zeros(Lx, Ly)
    mx_m .= mx
    my_m .= my

    Ms = 1e5
    phi = zeros(Nx, Nx)
    for i = 1:Nx, j=1:Ny
        phi[i, j] = NuMag.compute_magnetic_phase_direct(mx_m, my_m, dx, dy, Lz, Ms, i-80+16, j-80+32)
    end
    return phi
end

phi1 = phase_in_theory()
phi2 = phase_fft()
phi3 = phase_direct()
println(maximum(phi2.-phi1)/maximum(phi1))
println(maximum(phi3.-phi2)/maximum(phi2)< 1e-10)
@test maximum(phi2.-phi1)/maximum(phi1) < 3e-4