using JuMag
using Printf
using Test
using NPZ

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
    d = 1e-9
    a = F0(x-Lx/2, y-Ly/2) - F0(x+Lx/2, y-Ly/2) - F0(x-Lx/2, y+Ly/2) + F0(x+Lx/2, y+Ly/2)
    b = F0(y-Ly/2, x-Lx/2) - F0(y+Ly/2, x-Lx/2) - F0(y-Ly/2, x+Lx/2) + F0(y+Ly/2, x+Lx/2)
    return mu0*Ms*d^2*pi/(Phi0)*(my*b-mx*a)
end

function phase_in_theory(N)
    mx = cos(5/3*pi)
    my = sin(5/3*pi)
    Lx = 32
    Ly = 64
    Lz = 16
    Ms = 1e5
    phi = zeros(N, N)
    for i = 1:N, j=1:N
        x = (i - div(N, 2)) - 0.5
        y = (j - div(N, 2)) - 0.5
        phi[i, j] = phim_uniformly_magnetized_slab(x, y, mx, my, Lx, Ly, Lz, Ms)
    end
    return phi
end

function uniform_slab()
    nx, ny, nz = 32, 64, 16
    m = zeros(3,nx,ny,nz)
    m[1,:,:,:] .= cos(5/3*pi)
    m[2,:,:,:] .= sin(5/3*pi)
    return m
end

function phase_fft(N)
    m = uniform_slab()
    phi = JuMag.compute_magnetic_phase(m, 0.0, "X", N1=100, N2=N, Ms=1e5, d=1e-9)
    return phi
end

N=300
phi1 = phase_in_theory(N)
phi2 = phase_fft(N)
npzwrite("phi1.npy", phi1)
npzwrite("phi2.npy", phi2)

println(maximum(phi2.-phi1)/maximum(phi1))
@test maximum(phi2.-phi1)/maximum(phi1) < 3e-4

# using Plots
# N = 300
# ctof = 100
# phi1 = npzread("phi1.npy")[ctof:end-ctof,ctof:end-ctof]
# phi2 = npzread("phi2.npy")[ctof:end-ctof,ctof:end-ctof]
# p = plot(layout=grid(1,2))
# plot!(p[1],phi1,seriestype=:contour,ratio=:equal)
# plot!(p[2],phi2,seriestype=:contour,ratio=:equal)

# norm_phi1 = phi1/maximum(phi1)
# norm_phi2 = phi2/maximum(phi2)
# print(maximum(abs.(norm_phi1-norm_phi2)) / maximum(abs.(norm_phi1)))
