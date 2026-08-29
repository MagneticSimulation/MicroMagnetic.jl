using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using Random
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
    B = build_matrix(sim, gamma=2.21e5/mu_0, alpha=0.0)
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

function test_matrixfree_equivalence_afm(setup_fun; gamma=2.21e5/mu_0, alpha=0.0,
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

    rng = MersenneTwister(2024)
    N2 = size(B, 1)
    for _ in 1:nrand
        x = randn(rng, N2)
        y_ref = B * x
        y_mf  = similar(y_ref); LinearAlgebra.mul!(y_mf, op, x)
        @test isapprox(y_ref, y_mf; rtol=rtol, atol=atol)

        a5, b5 = 0.5, 3.1
        y_out = randn(rng, N2)
        y_5   = copy(y_out)
        LinearAlgebra.mul!(y_5, op, x, a5, b5)
        y_5_ref = a5 .* (B * x) .+ b5 .* y_out
        @test isapprox(y_5_ref, y_5; rtol=rtol, atol=atol)
    end
    nothing
end

H = 0.3
B, f = compute_frequency(H)
B_ex, fan = analytical(H)

println("f=", f, " ", fan)
@test abs(f -  fan)/f < 100*eps()
@test isapprox(-B, B_ex)

# --- matrix-free equivalence (2-site AFM + antiferro exch + uniaxial anis) --
println("  matrix-free: AFM (α=0) …")
test_matrixfree_equivalence_afm(() -> setup(H=(0, 0, H)))
println("  matrix-free: AFM damped (α=0.01) …")
test_matrixfree_equivalence_afm(() -> setup(H=(0, 0, H)); alpha=0.01)

