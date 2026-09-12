# Acceptance tests for the SweepMatch stage-chain scenario
# (src/adjoint/sweepmatch.jl).
#
# Protocol under test: rotation magnetization curve — large field
# (80 mT >> anisotropy field) rotating in the x–z plane through the easy axis,
# single stable branch per stage, 12 stages on a 12×4×1 grid.
#
#   (a) run_forward! reproduces a hand-written stage loop (same stage fields,
#       relax continuation, Ms-weighted observer, J accumulation);
#   (b) gradient! matches central finite differences for a Ku-map cell and a
#       scalar Aex (rel < 1e-4; stopping_dmdt tightened to 1e-9);
#   (c) path dependence: scrambling the checkpoint↔stage pairing changes the
#       gradient, and reversing a saturate/probe stage order swaps the
#       checkpointed states and hence J and the gradient;
#   (+) zero-weight intermediate stages contribute exactly nothing (their
#       ∇g projects to the zero vector, λ = 0, no solve) and the weighted
#       gradient equals the single-(last-)stage gradient;
#   (+) API contract: constructor/set_design!/gradient! guards.
using MicroMagnetic
using Test
using LinearAlgebra

# the shared test runner (test_utils.jl) leaves the last variant's global
# backend/precision in place (known mainline leak: e.g. fem/test_demag ends on
# CUDA); this file builds its own Float64 CPU sims, so pin both first.  The
# _n_sims zeroing only suppresses the "existing Sim(s)" warning.
_saved_n_sims = MicroMagnetic._n_sims[]
MicroMagnetic._n_sims[] = 0
MicroMagnetic.set_backend("cpu")
MicroMagnetic.set_precision(Float64)
MicroMagnetic._n_sims[] = _saved_n_sims

const NX, NY = 12, 4
const NSTAGE = 12
# exclude the duplicated endpoint of range(0, 2π; length=N+1): the protocol is
# one full turn without visiting phi = 0 twice
const ANGLES = collect(range(0, 2π; length = NSTAGE + 1)[1:(end - 1)])
const HMAG = 8e4                       # 80 mT >> anisotropy field (~16-25 mT)
const TGT = 0.9 .* cos.(ANGLES)        # off-curve targets => healthy losses
const KU_MAP = [5e3 + 3e3 * sin(π * i / (NX + 1)) for i in 1:NX for j in 1:NY]
const AEX0 = 1.3e-11
const CELL = 3 + (2 - 1) * NX          # x = 3, y = 2

function rotation_sim(; name, Ku, m0 = (0.6, 0.0, 0.8))
    mesh = FDMesh(nx = NX, ny = NY, nz = 1, dx = 5e-9, dy = 5e-9, dz = 5e-9)
    sim = Sim(mesh; driver = "SD", name = name)
    set_Ms(sim, 8e5)
    add_exch(sim, AEX0)
    add_anis(sim, Ku; axis = (0, 0, 1))
    init_m0(sim, m0; norm = true)
    return sim
end

function make_prob(; name = "sweepmatch", target = CurveTarget(Moment(:z), TGT))
    return SweepMatch(rotation_sim(name = name, Ku = KU_MAP);
                      design = (Ku = KU_MAP, Aex = AEX0),
                      sweep = (H_list = HMAG, angles = ANGLES),
                      target = target)
end

@testset "adjoint/sweepmatch" begin

    @testset "(a) forward reproduces a hand-written stage loop" begin
        prob = make_prob()
        J = run_forward!(prob)

        # hand loop: independent sim, same stage protocol, fresh SD phase per
        # stage (documented stage semantics), plain Ms-weighted observation
        sim = rotation_sim(name = "sweep_hand", Ku = KU_MAP)
        add_zeeman(sim, (HMAG * sin(ANGLES[1]), 0.0, HMAG * cos(ANGLES[1])); name = "hand_H")
        ms = Float64.(Array(sim.mu0_Ms))
        ms_sum = sum(ms)
        M_hand = Float64[]
        for phi in ANGLES
            update_zeeman(sim, (HMAG * sin(phi), 0.0, HMAG * cos(phi)); name = "hand_H")
            sim.driver.steps = 0                       # fresh Barzilai–Borwein phase
            relax(sim; max_steps = 20000, stopping_dmdt = 1e-9,
                  save_data_every = -1, save_m_every = -1)
            m = Float64.(Array(sim.spin))
            push!(M_hand, sum(ms[i] * m[3 * (i - 1) + 3] for i in 1:length(ms)) / ms_sum)
        end
        J_hand = sum((M_hand[k] - TGT[k])^2 for k in 1:NSTAGE)

        @test maximum(abs.(prob.obs_curve .- M_hand)) < 1e-7
        @test isapprox(J, J_hand; rtol = 1e-7, atol = 1e-12)
        # no relaxation-quality retry fired, so the two loops relaxed identically
        @test maximum(prob.relax_residual) < prob.relax_tol
        # checkpoint footprint: stages × (3N state + 3N field) Float64
        @test checkpoint_memory(prob) == NSTAGE * 2 * 3 * NX * NY * sizeof(Float64)
    end

    @testset "(b) gradient matches central finite differences" begin
        prob = make_prob(name = "sweep_fd")
        J0 = run_forward!(prob)
        G = gradient!(prob)

        # the adjoint solves themselves converged: relative residual < 1e-6,
        # except at degenerate stages whose ‖P_t∇g‖ sits at the zero short-circuit
        # (target ≈ observation to ~1e-13 ⇒ the relative residual amplifies
        # roundoff while the absolute error, bounded by ‖P_t∇g‖, is negligible)
        @test all(prob.adjoint_residual .< 1e-6 .|| prob.adjoint_bnorm .< 1e-10)
        # the Aex channel must be live: a uniform rotation would give exactly 0
        @test abs(G.Aex) > 1e3

        # Ku-map single cell
        eps = 1e-2 * KU_MAP[CELL]
        kp = copy(KU_MAP); kp[CELL] += eps
        set_design!(prob; Ku = kp); Jp = run_forward!(prob)
        km = copy(KU_MAP); km[CELL] -= eps
        set_design!(prob; Ku = km); Jm = run_forward!(prob)
        fd = (Jp - Jm) / (2eps)
        @test isapprox(fd, G.Ku[CELL]; rtol = 1e-4, atol = 1e-12)
        set_design!(prob; Ku = KU_MAP)

        # scalar Aex
        epsa = 1e-2 * AEX0
        set_design!(prob; Aex = AEX0 + epsa); Jp = run_forward!(prob)
        set_design!(prob; Aex = AEX0 - epsa); Jm = run_forward!(prob)
        fd = (Jp - Jm) / (2epsa)
        @test isapprox(fd, G.Aex; rtol = 1e-4, atol = 1e-12)
        set_design!(prob; Aex = AEX0)

        # after restoring the design the forward is reproducible
        @test run_forward!(prob) == J0
    end

    @testset "(+) zero-weight stages: λ = projection of the zero vector" begin
        w = zeros(NSTAGE); w[end] = 1.0
        prob = make_prob(name = "sweep_w",
                         target = CurveTarget(Moment(:z), TGT; weights = w))
        J = run_forward!(prob)
        G = gradient!(prob)

        @test isapprox(J, (prob.obs_curve[end] - TGT[end])^2; rtol = 1e-12)
        # zero-weight stages short-circuit: no solve, λ = 0, residual 0
        @test all(prob.adjoint_residual[1:(end - 1)] .== 0.0)
        @test prob.adjoint_residual[end] < 1e-6

        # equals the gradient of the single-(last-)stage problem: both stage
        # states are the same unique equilibrium at Ĥ(end)
        prob1 = SweepMatch(rotation_sim(name = "sweep_1stage", Ku = KU_MAP);
                           design = (Ku = KU_MAP, Aex = AEX0),
                           sweep = (H_list = HMAG, angles = [ANGLES[end]]),
                           target = CurveTarget(Moment(:z), [TGT[end]]))
        run_forward!(prob1)
        G1 = gradient!(prob1)
        @test norm(G.Ku - G1.Ku) / norm(G.Ku) < 1e-2
        @test isapprox(G.Aex, G1.Aex; rtol = 1e-2)
    end

    @testset "(c) path dependence / checkpoint order" begin
        prob = make_prob(name = "sweep_path")
        run_forward!(prob)
        G = gradient!(prob)

        # (i) scramble the checkpoint↔stage pairing: stage k is solved at the
        # state of stage k+1 — a different (wrong) inverse problem
        prob.cp_m = circshift(prob.cp_m, -1)
        prob.cp_H0 = circshift(prob.cp_H0, -1)
        Gs = gradient!(prob)
        @test norm(Gs.Ku .- G.Ku) / norm(G.Ku) > 0.1
        @test norm(Gs.Aex .- G.Aex) / norm(G.Aex) > 0.1

        # (ii) reverse a saturate/probe stage order: the probe field (half the
        # anisotropy field, opposing, slightly tilted) holds the moment in the
        # +z basin on a canted metastable branch while the saturate stage pins
        # it to exactly +z — so the two paths visit the SAME two states in
        # SWAPPED stage slots, which moves J (asymmetric targets) and the
        # gradient
        tilt = deg2rad(10)
        d_sat = (0.0, 0.0, 1.0)                # saturate exactly along the axis
        d_probe = (sin(tilt), 0.0, -cos(tilt)) # opposing, slightly tilted
        hk = 2 * 5e4 / (4e-7 * pi * 8e5)       # anisotropy field, A/m
        probe = 0.5 * hk                       # sub-switching opposing magnitude
        simA = rotation_sim(name = "sweep_pa", Ku = 5e4, m0 = (0.05, 0.0, 1.0))
        simB = rotation_sim(name = "sweep_pb", Ku = 5e4, m0 = (0.05, 0.0, 1.0))
        pa = SweepMatch(simA; design = (Ku = 5e4,),
                        sweep = (H_list = [8e5, probe],
                                 directions = [d_sat, d_probe]),
                        target = CurveTarget(Moment(:z), [0.5, 0.95]))
        pb = SweepMatch(simB; design = (Ku = 5e4,),
                        sweep = (H_list = [probe, 8e5],
                                 directions = [d_probe, d_sat]),
                        target = CurveTarget(Moment(:z), [0.5, 0.95]))
        Ja = run_forward!(pa); Ga = gradient!(pa)
        Jb = run_forward!(pb); Gb = gradient!(pb)

        @test pa.obs_curve[1] ≈ pb.obs_curve[2] atol = 1e-6    # saturate stages
        @test pa.obs_curve[2] ≈ pb.obs_curve[1] atol = 1e-6    # probe stages
        @test pa.obs_curve[1] - pb.obs_curve[1] > 0.005        # genuinely swapped
        @test abs(Ja - Jb) > 1e-3
        @test abs(Ga.Ku - Gb.Ku) / max(abs(Ga.Ku), abs(Gb.Ku)) > 0.05
    end

    @testset "API contract" begin
        sim = rotation_sim(name = "sweep_api", Ku = KU_MAP)
        @test_throws ArgumentError SweepMatch(sim; design = (Ms = 8e5,),
                                              sweep = (H_list = HMAG, angles = ANGLES),
                                              target = CurveTarget(Moment(:z), TGT))
        @test_throws ArgumentError Moment(:w)
        @test_throws ArgumentError SweepMatch(sim; design = (Ku = KU_MAP,),
                                              sweep = (H_list = HMAG,),
                                              target = CurveTarget(Moment(:z), TGT))
        @test_throws ArgumentError SweepMatch(sim; design = (Ku = KU_MAP,),
                                              sweep = (H_list = HMAG, angles = ANGLES,
                                                       directions = [(0, 0, 1)]),
                                              target = CurveTarget(Moment(:z), TGT))
        @test_throws DimensionMismatch SweepMatch(sim; design = (Ku = KU_MAP,),
                                                  sweep = (H_list = HMAG, angles = ANGLES),
                                                  target = CurveTarget(Moment(:z),
                                                                       fill(0.5, NSTAGE + 1)))
        prob = SweepMatch(sim; design = (Ku = KU_MAP,),
                          sweep = (H_list = HMAG, angles = ANGLES),
                          target = CurveTarget(Moment(:z), TGT))
        @test_throws ErrorException gradient!(prob)          # no forward yet
        run_forward!(prob)
        @test gradient!(prob) isa NamedTuple
        set_design!(prob; Ku = KU_MAP .* 2)                   # invalidates forward
        @test_throws ErrorException gradient!(prob)
        @test_throws ArgumentError set_design!(prob; Aex = 1e-11)   # undeclared
        @test_throws ArgumentError set_design!(prob; Ku = 5e3)      # map got scalar
        @test_throws ArgumentError SweepMatch(sim; design = (Ku = KU_MAP,),
                                              sweep = (H_list = HMAG, angles = ANGLES),
                                              target = CurveTarget(Moment(:z), TGT),
                                              zeeman_name = "sweep_H")  # duplicate
    end
end
