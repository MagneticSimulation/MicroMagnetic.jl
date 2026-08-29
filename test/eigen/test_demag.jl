using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using Random
using Test

MicroMagnetic.set_backend("cpu")

function setup(;m0=(1,0,0), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=1e-9, dz=1e-9)
    
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, m0)
    add_zeeman(sim, H)
    add_demag(sim)
    return sim
end

MicroMagnetic.set_precision(AbstractFloat)

function compute_frequency(H0)
    H = (H0, 0, 0)
    sim = setup(H=H)
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
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

function test_matrixfree_equivalence_demag(setup_fun; gamma=2.21e5, alpha=0.0,
                                            rtol=1e-12, atol=1e-14, nrand=5)
    _α = Float64(alpha)
    sim_d = setup_fun()
    sim_o = setup_fun()
    B  = build_matrix(sim_d, gamma=Float64(gamma), alpha=_α)
    op = build_matrix(sim_o, gamma=Float64(gamma), alpha=_α, matrixfree=true)
    @assert op.alpha ≈ _α "matrix-free operator has alpha=$(op.alpha), want $_α"

    @test size(op) == size(B)
    @test eltype(op) == eltype(B)

    B_mf = Matrix(op)
    @test isapprox(B, B_mf; rtol=rtol, atol=atol)

    rng = MersenneTwister(9999)
    N2 = size(B, 1)
    for _ in 1:nrand
        x = randn(rng, N2)
        y_ref = B * x
        y_mf  = similar(y_ref); LinearAlgebra.mul!(y_mf, op, x)
        @test isapprox(y_ref, y_mf; rtol=rtol, atol=atol)

        a5, b5 = -2.0, 0.25
        y_out = randn(rng, N2)
        y_5   = copy(y_out)
        LinearAlgebra.mul!(y_5, op, x, a5, b5)
        y_5_ref = a5 .* (B * x) .+ b5 .* y_out
        @test isapprox(y_5_ref, y_5; rtol=rtol, atol=atol)
    end
    nothing
end

H = 1.23e4
B, f = compute_frequency(H)
fan = analytical(H)

println("f=", f, " ", fan)
@test abs(f -  fan)/f < 100*eps()

# --- matrix-free equivalence (single-cell demag + Zeeman, long-range path) ---
println("  matrix-free: 1-cell demag (α=0) …")
test_matrixfree_equivalence_demag(() -> setup(H=(H, 0, 0)))
println("  matrix-free: 1-cell demag damped (α=0.01) …")
test_matrixfree_equivalence_demag(() -> setup(H=(H, 0, 0)); alpha=0.01)


