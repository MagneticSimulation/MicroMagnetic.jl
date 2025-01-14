using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using Test

MicroMagnetic.set_backend("cpu")

function setup(;m0=(0,0,1), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=2e-9)
    sim = create_sim(mesh; H=H, m0=m0, Ms=8e5, Kc=2e4, demag=false)
    return sim
end

MicroMagnetic.set_precision(AbstractFloat)

function compute_frequency100(H0)
    H = (H0, 0, 0)
    sim = setup(m0=(1,0,0), H=H)
    B = build_matrix(sim, gamma=2.21e5)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

function compute_frequency110(H0)
    H = (H0/sqrt(2), H0/sqrt(2),0)
    sim = setup(H=H, m0=(1,1,0))
    B = build_matrix(sim, gamma=2.21e5)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

function analytical100(H, Kc=2e4, Ms=8e5)
    gamma = 2.21e5
    K = 4*Kc/(mu_0*Ms)
    freq = gamma*(H .+ K)/1e9/(2*pi)
    return freq
end

function analytical110(H, Kc=2e4, Ms=8e5)
    gamma = 2.21e5
    K = Kc/(mu_0*Ms)
    freq = gamma*sqrt.(H.^2 .- 2K.*H .-8K^2)/1e9/(2*pi)
    return freq
end

H = 1e5
f100 = compute_frequency100(H)
f100_an = analytical100(H)

f110 = compute_frequency110(H)
f110_an = analytical110(H)

println("f100:", f100, " ", f100_an)
println("f110:", f110, " ", f110_an)
@test f100 == f100_an
@test f110 == f110_an

using Enzyme
MicroMagnetic.set_precision(Float64)

function compute_frequency100_enzyme(H0)
    H = (H0, 0, 0)
    sim = setup(m0=(1,0,0), H=H)
    B = dynamic_matrix(sim, gamma=2.21e5)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

function compute_frequency110_enzyme(H0)
    H = (H0/sqrt(2), H0/sqrt(2),0)
    sim = setup(H=H, m0=(1,1,0))
    B = dynamic_matrix(sim, gamma=2.21e5)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

f100 = compute_frequency100_enzyme(H)
f110 = compute_frequency110_enzyme(H)
println("f100_enzyme:", f100, " ", f100_an)
println("f110_enzyme:", f110, " ", f110_an)
@test abs(f100-f100_an) < eps()
@test abs(f110-f110_an) < eps()*10
