# Acceptance tests for the FrequencyMatch scenario.
#
# Double acceptance:
#   (a) analytical: uniform uniaxial FMR ω = γ(H + 2Ku/(μ0·Ms)) — the
#       selected mode, and the chained dJ/dKu from gradient!, vs the closed form
#       (relative error < 1e-6; GHz convention as in test/eigen/test_cubic.jl);
#   (b) finite differences: perturb the design variable, rebuild B, re-solve —
#       central-difference dω/dθ vs analytical vs adjoint, three-way cross-check
#       (Ku on a single cell; Aex and D on a 4-cell exchange/DMI chain where the
#       selected k=1 mode makes both derivatives non-trivial).

using MicroMagnetic
using Test
using LinearAlgebra
using Random

MicroMagnetic.set_backend("cpu")
MicroMagnetic.set_precision(AbstractFloat)

const GAMMA = 2.21e5
freq_of(lambda) = imag(lambda) / (2π * 1e9)     # GHz, test_cubic convention

function cell_sim(; Ku, H, Ms = 8e5)
    mesh = FDMesh(nx = 1, ny = 1, nz = 1, dx = 5e-9, dy = 5e-9, dz = 2e-9)
    sim = Sim(mesh; driver = "SD")
    set_Ms(sim, Ms)
    add_anis(sim, Ku)
    H == 0 || add_zeeman(sim, (0, 0, H))
    init_m0(sim, (0, 0, 1))
    return sim
end

function chain_sim(; A, Ku, D = 0.0, axis = (0, 0, 1))
    mesh = FDMesh(nx = 4, ny = 1, nz = 1, dx = 5e-9, dy = 5e-9, dz = 5e-9)
    sim = Sim(mesh; driver = "SD")
    set_Ms(sim, 8e5)
    add_exch(sim, A)
    add_anis(sim, Ku; axis = axis)
    D == 0 || add_dmi(sim, D; type = "bulk")
    init_m0(sim, axis)
    return sim
end

# uniform uniaxial FMR: ω = γ(H + 2Ku/(μ0·Ms)) [rad/s] → GHz
omega_ana(Ku; H = 8e4, Ms = 8e5, gamma = GAMMA) =
    gamma * (H + 2 * Ku / (mu_0 * Ms)) / (2π * 1e9)
domega_dKu_ana(; Ms = 8e5, gamma = GAMMA) = gamma * 2 / (mu_0 * Ms) / (2π * 1e9)

@testset "adjoint/freqmatch" begin

    @testset "analytical: uniform uniaxial FMR" begin
        Ms, Ku0, H0 = 8e5, 5e4, 8e4
        ω0 = omega_ana(Ku0; H = H0, Ms = Ms)
        ωstar = ω0 + 0.05          # off-resonance target → non-trivial chain rule

        sim = cell_sim(Ku = Ku0, H = H0, Ms = Ms)
        fm = FrequencyMatch(sim; design = (Ku = Ku0,), target = (freq = ωstar),
                            gamma = GAMMA, alpha = 0.0)
        J = run_forward!(fm)

        @test fm.omega > 0                      # ±iω tie resolved to the +Im branch
        @test fm.omega ≈ ω0 rtol = 1e-9         # selected mode = analytical FMR
        @test J ≈ (fm.omega - ωstar)^2 rtol = 1e-12

        G = gradient!(fm)
        dJ_dKu_ana = 2 * (fm.omega - ωstar) * domega_dKu_ana(Ms = Ms)
        @test G.Ku ≈ dJ_dKu_ana rtol = 1e-6
        @test keys(G) == (:Ku,)
    end

    @testset "FD cross-check: Ku (analytical vs FD vs adjoint)" begin
        Ms, Ku0, H0 = 8e5, 5e4, 8e4
        ωstar = omega_ana(Ku0; H = H0, Ms = Ms) + 0.1

        sim = cell_sim(Ku = Ku0, H = H0, Ms = Ms)
        fm = FrequencyMatch(sim; design = (Ku = Ku0,), target = (freq = ωstar),
                            gamma = GAMMA, alpha = 0.0)
        run_forward!(fm)
        G = gradient!(fm)
        dω_adj = G.Ku / (2 * (fm.omega - ωstar))
        dω_ana = domega_dKu_ana(Ms = Ms)

        δ = 1.0                                  # ω is linear in Ku → FD exact to roundoff
        set_design!(fm; Ku = Ku0 + δ); run_forward!(fm); ωp = fm.omega
        set_design!(fm; Ku = Ku0 - δ); run_forward!(fm); ωm = fm.omega
        dω_fd = (ωp - ωm) / (2δ)

        @test dω_adj ≈ dω_ana rtol = 1e-6
        @test dω_fd ≈ dω_ana rtol = 1e-6
        @test dω_fd ≈ dω_adj rtol = 1e-6

        set_design!(fm; Ku = Ku0); run_forward!(fm)   # round trip restores the point
        @test fm.omega ≈ omega_ana(Ku0; H = H0, Ms = Ms) rtol = 1e-9
    end

    @testset "FD cross-check: Aex (4-cell chain, k=1 mode)" begin
        A0, Ku0 = 1.3e-11, 5e4
        sim = chain_sim(A = A0, Ku = Ku0)
        # chain modes: k=0 at ~3.5 GHz (exchange-blind), k=1 at ~25 GHz, k=2 at ~76 GHz
        fm = FrequencyMatch(sim; design = (Aex = A0, Ku = Ku0), target = (freq = 25.0),
                            gamma = GAMMA, alpha = 0.0)
        run_forward!(fm)
        @test 20 < fm.omega < 30                 # the k=1 mode, not k=0
        ω0sel, ωstar = fm.omega, 25.0

        G = gradient!(fm)
        dω_adj = G.Aex / (2 * (ω0sel - ωstar))

        δ = 1e-12                                # ω linear in Aex → FD exact to roundoff
        set_design!(fm; Aex = A0 + δ); run_forward!(fm); ωp = fm.omega
        set_design!(fm; Aex = A0 - δ); run_forward!(fm); ωm = fm.omega
        dω_fd = (ωp - ωm) / (2δ)

        @test abs(dω_adj) > 1e-3                 # exchange shifts the k=1 mode
        @test dω_fd ≈ dω_adj rtol = 1e-6
    end

    @testset "FD cross-check: bulk D (4-cell in-plane chain)" begin
        A0, Ku0, D0 = 1.3e-11, 5e4, 2e-4
        # in-plane ground state (easy axis x̂): only here does the bulk-DMI chiral
        # term enter the linearisation at all — with m0 = ẑ the tangent plane
        # kills δm_z to first order and ∂B/∂D vanishes identically.  The chiral
        # shift is even in D, so the derivative is taken at D0 ≠ 0.
        sim = chain_sim(A = A0, Ku = Ku0, D = D0, axis = (1, 0, 0))
        fm = FrequencyMatch(sim; design = (D = D0,), target = (freq = 27.5),
                            gamma = GAMMA, alpha = 0.0)
        run_forward!(fm)
        @test 20 < fm.omega < 30                 # the k=1 branch, away from k=0 / k=2
        ω0sel, ωstar = fm.omega, 27.5

        G = gradient!(fm)
        dω_adj = G.D / (2 * (ω0sel - ωstar))

        δ = 1e-6
        set_design!(fm; D = D0 + δ); run_forward!(fm); ωp = fm.omega
        set_design!(fm; D = D0 - δ); run_forward!(fm); ωm = fm.omega
        dω_fd = (ωp - ωm) / (2δ)

        @test abs(dω_adj) > 1e-3                 # DMI moves the k=1 branch
        @test dω_fd ≈ dω_adj rtol = 1e-6
    end

    @testset "Adjoint(op) mul! agrees with dense Bᵀ" begin
        sim = cell_sim(Ku = 5e4, H = 8e4)
        fm = FrequencyMatch(sim; design = (Ku = 5e4,), target = (freq = 10.0),
                            gamma = GAMMA, alpha = 0.0)
        run_forward!(fm)

        # left/right eigenvector residuals (dense B from the forward cache)
        B, λ, x, y = fm.B, fm.lambda, fm.xsel, fm.ysel
        scale = 1.0 + abs(λ)
        @test norm(B * x - λ * x) <= 1e-8 * scale
        @test norm(B' * y - λ * y) <= 1e-8 * scale

        # the P1 matrix-free adjoint reproduces Bᵀ on random vectors — the same
        # operator the left eigenvector machinery will ride on for P5/Arpack
        op = build_matrix(fm.sim; gamma = GAMMA, alpha = 0.0, matrixfree = true)
        rng = MersenneTwister(2026)
        for _ in 1:3
            v = randn(rng, size(B, 1))
            out = zeros(size(B, 1))
            mul!(out, adjoint(op), v)
            @test out ≈ B' * v rtol = 1e-10
        end
    end

    @testset "protocol: set_design!, value_and_gradient, auto-forward" begin
        ωstar = omega_ana(5e4) + 0.05
        sim = cell_sim(Ku = 5e4, H = 8e4)
        fm = FrequencyMatch(sim; design = (Ku = 5e4,), target = ωstar,   # bare-Real target
                            gamma = GAMMA, alpha = 0.0)

        G = gradient!(fm)                        # no explicit run_forward! first
        @test fm.J !== nothing
        @test G.Ku ≈ 2 * (fm.omega - ωstar) * domega_dKu_ana() rtol = 1e-6

        J1 = run_forward!(fm)
        set_design!(fm; Ku = 6e4)
        @test fm.J === nothing                   # cache invalidated
        @test fm.design.Ku == 6e4
        J2 = run_forward!(fm)
        @test J2 ≈ (omega_ana(6e4) - ωstar)^2 rtol = 1e-9
        @test J2 != J1

        Jv, Gv = value_and_gradient(fm)
        @test Jv ≈ (omega_ana(6e4) - ωstar)^2 rtol = 1e-9
        @test Gv isa NamedTuple && keys(Gv) == (:Ku,)
    end

    @testset "error paths" begin
        sim = cell_sim(Ku = 5e4, H = 8e4)
        @test_throws ArgumentError FrequencyMatch(sim; design = (Ms = 8e5,),
                                                  target = (freq = 1.0))
        @test_throws ArgumentError FrequencyMatch(sim; design = (Aex = 1.3e-11,),
                                                  target = (freq = 1.0))   # no exchange in sim
        @test_throws ArgumentError FrequencyMatch(sim; design = (Ku = 5e4, Foo = 1.0),
                                                  target = (freq = 1.0))
        @test_throws ArgumentError FrequencyMatch(sim; design = (Ku = [5e4, 6e4],),
                                                  target = (freq = 1.0))
        @test_throws ArgumentError FrequencyMatch(sim; design = (Ku = 5e4,),
                                                  target = :bad)             # not a frequency
        @test_throws ArgumentError set_design!(FrequencyMatch(sim; design = (Ku = 5e4,),
                                               target = (freq = 1.0), gamma = GAMMA); D = 1e-3)

        mesh = FDMesh(nx = 1, ny = 1, nz = 1)    # no anisotropy → :Ku unresolvable
        sim2 = Sim(mesh; driver = "SD")
        set_Ms(sim2, 8e5)
        init_m0(sim2, (0, 0, 1))
        @test_throws ArgumentError FrequencyMatch(sim2; design = (Ku = 5e4,),
                                                  target = (freq = 1.0))
    end
end
