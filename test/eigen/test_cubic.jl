using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using Random
using Test

MicroMagnetic.set_backend("cpu")
MicroMagnetic.set_precision(AbstractFloat)

function setup(;m0=(0,0,1), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=2e-9)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5)
    add_cubic_anis(sim, 2e4)
    add_zeeman(sim, H)
    init_m0(sim, m0)
    return sim
end

function compute_frequency100(H0)
    H = (H0, 0, 0)
    sim = setup(m0=(1,0,0), H=H)
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

function compute_frequency110(H0)
    H = (H0/sqrt(2), H0/sqrt(2),0)
    sim = setup(H=H, m0=(1,1,0))
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
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

# Verify matrix-free operator B_mf ≡ dense B point-by-point and via mul! on
# random tangent vectors (both 3- and 5-argument forms required by Arpack).
function test_matrixfree_equivalence(setup_fun; gamma=2.21e5, alpha=0.0,
                                     rtol=1e-12, atol=1e-14, nrand=5)
    _α = Float64(alpha)
    sim_d = setup_fun()
    sim_o = setup_fun()
    # Comma-only kwargs (no semicolon) matches the compute_frequency* call
    # style that has been verified to pass the analytical checks above.
    B  = build_matrix(sim_d, gamma=Float64(gamma), alpha=_α)
    op = build_matrix(sim_o, gamma=Float64(gamma), alpha=_α, matrixfree=true)

    # Guard against keyword mis-dispatch (the alpha=0.01 default must NOT
    # kick in when the caller requests alpha=0.0).
    @assert op.alpha ≈ _α "matrix-free operator has alpha=$(op.alpha), want $_α"
    if size(B, 1) == 2
        # For 1-site systems the damping ratio α = |λ/ω| from eigenvalues
        # gives a direct fingerprint of the effective α used by the dense path.
        v = eigvals(B)
        ω, λ = maximum(abs.(imag(v))), maximum(abs.(real(v)))
        if ω > 0 && abs(_α) < eps()
            @assert iszero(λ) || λ / ω < 1e-4 "dense build: damping λ/ω=$(λ/ω) ≠ 0 for alpha=0 (possible keyword-mixup → alpha defaulted to 0.01)"
        end
    end

    @test size(op) == size(B)
    @test eltype(op) == eltype(B)
    @test issymmetric(op) == false
    @test ishermitian(op) == false

    # Full-column materialisation: the two Jacobians must match exactly.
    B_mf = Matrix(op)
    @test isapprox(B, B_mf; rtol=rtol, atol=atol)

    # 3-arg mul!: op * x vs B * x on random tangent vectors.
    rng = MersenneTwister(1234)
    N2 = size(B, 1)
    for _ in 1:nrand
        x = randn(rng, N2)
        y_ref = B * x
        y_mf  = similar(y_ref); LinearAlgebra.mul!(y_mf, op, x)
        @test isapprox(y_ref, y_mf; rtol=rtol, atol=atol)

        # 5-arg mul!: y = α*B*x + β*y (Arpack hot path)
        a5 = 2.3
        b5 = -1.7
        y_out = randn(rng, N2)
        y_5   = copy(y_out)
        LinearAlgebra.mul!(y_5, op, x, a5, b5)
        y_5_ref = a5 .* (B * x) .+ b5 .* y_out
        @test isapprox(y_5_ref, y_5; rtol=rtol, atol=atol)
    end

    # copy() must produce an independent op (different caches, same result).
    op2 = copy(op)
    y2 = zeros(N2); x2 = randn(rng, N2)
    LinearAlgebra.mul!(y2, op2, x2)
    @test isapprox(y2, B * x2; rtol=rtol, atol=atol)
    @test op2._dm !== op._dm
    nothing
end

H = 1e5
f100 = compute_frequency100(H)
f100_an = analytical100(H)

f110 = compute_frequency110(H)
f110_an = analytical110(H)

println("f100:", f100, " ", f100_an)
println("f110:", f110, " ", f110_an)
# Eigenvalue extraction + symbolic → F64 roundtrip can differ in the last
# 1–2 ULPs from the analytic expression; a strict `==` is therefore flaky.
# Use a tight relative tolerance instead.
@test isapprox(f100, f100_an; rtol=10 * eps())
@test isapprox(f110, f110_an; rtol=10 * eps())

# --- matrix-free equivalence tests (cubic anisotropy + Zeeman, local-only) --
println("  matrix-free: [100] direction …")
test_matrixfree_equivalence(() -> setup(m0=(1,0,0), H=(H, 0, 0)))
println("  matrix-free: [110] direction …")
test_matrixfree_equivalence(() -> setup(m0=(1,1,0),
                                        H=(H/sqrt(2), H/sqrt(2), 0)))
# With non-zero damping: exercises the α·(m0×H0) baseline path too.
println("  matrix-free: damped α=0.01 …")
test_matrixfree_equivalence(() -> setup(m0=(1,0,0), H=(H, 0, 0)); alpha=0.01)