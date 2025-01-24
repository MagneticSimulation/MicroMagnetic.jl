using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using Test

MicroMagnetic.set_backend("cpu")

function setup(;m0=(1,0,0), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=1e-9, dz=1e-9)
    
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, m0)
    add_zeeman(sim, H)
    add_demag(sim, fft=false)
    return sim
end

MicroMagnetic.set_precision(AbstractFloat)

function compute_frequency(H0)
    H = (H0, 0, 0)
    sim = setup(H=H)
    B = build_matrix(sim, gamma=2.21e5)
    return B, imag(eigvals(B)[2])/1e9/(2*pi)
end

function analytical(H, Ms=8e5)
    Nx = 0.08831574004542228
    Ny = (1-Nx)/2
    K = 1/2*(Ny-Nx)*mu_0*Ms^2
    gamma = 2.21e5
    K = 2*K/(mu_0*Ms)
    freq = gamma*(H .+ K)/1e9/(2*pi)
    return freq
end

H = 1.23e4
B, f = compute_frequency(H)
fan = analytical(H)

println("f=", f, " ", fan)
@test abs(f -  fan)/f < 100*eps()


using Enzyme
MicroMagnetic.set_precision(Float64)

function compute_frequency_enzyme(H0)
    H = (H0, 0, 0)
    sim = setup(H=H)
    
    B = dynamic_matrix(sim, gamma=2.21e5)
    return B, imag(eigvals(B)[2])/1e9/(2*pi)
end

B2, fen = compute_frequency_enzyme(H)
println("fen=", fen, "  ", fan)

@test abs(fen -  fan)/f < 100*eps()

@test isapprox(B2, -B)

