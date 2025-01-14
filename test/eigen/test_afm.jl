using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using Test

MicroMagnetic.set_backend("cpu")

function spatial_m0(i, j, k, dx, dy, dz)
    if i == 1
        return (0,0,1)
    else
        return (0,0,-1)
    end
end

function setup(;H=(0,0,0))
    mesh = CubicMesh(nx=2, ny=1, nz=1)

    sim = Sim(mesh)
    set_mu_s(sim, mu_B)

    init_m0(sim, spatial_m0)

    J = 1meV
    add_exch(sim, -J)
    add_anis(sim, 0.01*J; axis=(0, 0, 1))
    add_zeeman(sim, H)
    return sim
end

MicroMagnetic.set_precision(AbstractFloat)

function compute_frequency(H0)
    H = (0, 0, H0)
    sim = setup(H=H)
    B = build_matrix(sim, gamma=2.21e5/mu_0)
    return B, imag(eigvals(B)[2])/1e9/(2*pi)
end


function analytical(H)
    gamma = 2.21e5/mu_0
    J = 1meV
    K = 0.01*J
    mu_s = mu_B
    we = gamma*J/mu_s
    wa = 2*gamma*K/mu_s
    wh = gamma*H
    w1 = wa+we+wh
    w2 = wa+we-wh
    A = [0  w1  0 we; -w1 0 we 0; 0 we 0 w2; we 0 -w2 0]
    freq = wh+sqrt(wa^2+2we*wa)
    return A, freq/1e9/(2*pi)
end

H = 0.3
B, f = compute_frequency(H)
B_ex, fan = analytical(H)

println("f=", f, " ", fan)
@test abs(f -  fan)/f < 100*eps()
@test isapprox(-B, B_ex)


using Enzyme
MicroMagnetic.set_precision(Float64)

function compute_frequency_enzyme(H0)
    H = (0, 0, H0)
    sim = setup(H=H)
    
    B = dynamic_matrix(sim, gamma=2.21e5/mu_0)
    return B, imag(eigvals(B)[2])/1e9/(2*pi)
end

B2, fen = compute_frequency_enzyme(H)
println("fen=", fen, "  ", fan)

@test isapprox(B2, B_ex)

@test abs(fen -  fan)/f < 100*eps()
