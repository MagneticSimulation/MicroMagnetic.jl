# StateMatch end-to-end FD validation.
#
# Central finite differences of J(θ) — each J a full forward (fixed init_m0 +
# relax to dmdt ≤ 1e-9) — against gradient!'s adjoint, per design parameter:
# scalar Ku / Aex / D(bulk), a single-cell Ku-map perturbation, a region
# parameterization (set_region + region_map, exercising the exchange/DMI pair
# table invalidation on write-back), and a vacuum-strip mesh (masked loss,
# vacuum-tolerant adjoint solve).
#
# Forward-convergence notes baked into the setup below (see StateMatch
# docstring): precession off (same equilibria, no fast-precession step pinning)
# and RK tolerance 1e-13 so the dmdt=1e-9 stopping floor is reachable.
using Test
using MicroMagnetic
using Printf

MicroMagnetic.set_backend("cpu")

mktempdir() do tmp
    cd(tmp) do
        nx = ny = 16
        Ms0 = 8e5
        Aex0 = 1.3e-11
        Ku0 = 5e4
        D0 = 1.5e-3
        alpha = 1.0

        # infeasible uniform-tilt target: J > 0 and every design gradient live
        # (the bulk DMI texture keeps m* non-uniform).  D0 stays at 1.5e-3
        # (D·d/A ≈ 0.58): near D ≈ 3e-3 the relaxed texture switches branches
        # across the FD interval (bifurcation), which poisons central
        # differences — the adjoint then equals the one-sided branch gradient.
        target(i, j, k, dx, dy, dz) = (sin(0.5), 0.0, cos(0.5))

        function base_sim(; ku = Ku0, dmi = true, ms = Ms0)
            mesh = FDMesh(nx = nx, ny = ny, nz = 1, dx = 5e-9, dy = 5e-9, dz = 5e-9)
            sim = Sim(mesh)
            set_Ms(sim, ms)
            add_exch(sim, Aex0)
            add_anis(sim, ku; axis = (0, 0, 1))
            dmi && add_dmi(sim, D0; type = "bulk")
            add_zeeman(sim, (1e5, 0, 0))
            set_alpha(sim, alpha)
            sim.driver.precession = false
            sim.driver.integrator.tol = 1e-13
            return sim
        end

        function make_sm(sim; design)
            StateMatch(sim; design = design, target_m = target, init_m = (0, 0, 1),
                       max_steps = 150000, stopping_dmdt = 1e-9, tol = 1e-10)
        end

        # FD harness: adjoint gradient at the base design, then J(θ±δ) forwards.
        # `shift!` moves the current design value RELATIVELY (through
        # set_design! — the production write-back path, pair-table drops
        # included); `extract` reads the G component.
        function fd_check(name, sm, extract, shift!, δ; rtol = 1e-5)
            run_forward!(sm)
            G = gradient!(sm)
            g_ad = Float64(extract(G))

            shift!(+δ);  Jp = run_forward!(sm)
            shift!(-2δ); Jm = run_forward!(sm)
            shift!(+δ)                            # restore
            g_fd = (Jp - Jm) / (2δ)

            err = abs(g_ad - g_fd) / max(abs(g_fd), abs(g_ad), 1e-30)
            println(rpad("  ✓ " * name, 34), "  adjoint = ", @sprintf("%+.6e", g_ad),
                    "   fd = ", @sprintf("%+.6e", g_fd),
                    "   rel err = ", @sprintf("%.2e", err))
            @test err < rtol
            return g_ad
        end

        @testset "adjoint/statematch FD" begin

            # ---- shared sim: scalar Ku + Aex + D designs (DMI-textured m*) ----
            sim = base_sim()
            sm = make_sm(sim; design = (Ku = Ku0, Aex = Aex0, D = D0))

            @testset "scalar designs" begin
                ku = Ref(Ku0)
                fd_check("scalar Ku", sm, G -> G.Ku, s -> (ku[] += s;
                                                           set_design!(sm; Ku = ku[])), 50.0)
                aex = Ref(Aex0)
                fd_check("scalar Aex", sm, G -> G.Aex, s -> (aex[] += s;
                                                             set_design!(sm; Aex = aex[])),
                         2.6e-16)
                d = Ref(D0)
                fd_check("scalar D (bulk)", sm, G -> G.D, s -> (d[] += s;
                                                                set_design!(sm; D = d[])),
                         3e-6)
            end

            # both solvers agree on the same frozen forward (DMI case)
            @testset "backend agreement" begin
                run_forward!(sm)
                g_kr = gradient!(sm).Ku
                # :pseudotime is stiffness-limited on textured states (the
                # adaptive step is pinned by the fast modes), so give it a
                # looser tolerance and more iterations for the comparison
                sm.adjoint_backend = :pseudotime
                sm.adjoint_tol = 1e-8
                sm.adjoint_maxiter = 30000
                g_pt = gradient!(sm).Ku
                sm.adjoint_backend = :krylov
                sm.adjoint_tol = 1e-10
                err = abs(g_kr - g_pt) / max(abs(g_kr), 1e-30)
                println("  ✓ krylov vs pseudotime (G.Ku)          rel diff = ",
                        @sprintf("%.2e", err))
                @test err < 1e-4
            end

            # ---- Ku map: single-cell perturbation (inline stencil path) ----
            @testset "Ku map (single cell)" begin
                ku_map = [Ku0 + 1.5e4 * sin(0.7 * i + 0.3 * j) for i in 1:nx
                          for j in 1:ny]
                sim2 = base_sim()
                sm2 = make_sm(sim2; design = (Ku = ku_map,))
                c0 = 8 + (8 - 1) * nx          # interior cell (i=8, j=8)
                k = copy(ku_map)
                fd_check("Ku map cell (8,8)", sm2, G -> G.Ku[c0],
                         s -> (k[c0] += s; set_design!(sm2; Ku = k)),
                         1e3)
            end

            # ---- region parameterization: 2 half-plane regions, region_map ----
            # design value is a function → per-cell gradient; the region
            # aggregate is the sum over the region's cells.  Region-constant Ku
            # takes the stencil pair-table fast path, so this also validates
            # set_design!'s cache invalidation (stale tables would corrupt J±).
            @testset "region parameterization" begin
                sim3 = base_sim()
                mesh = sim3.mesh                      # left / right half planes
                # (FDMesh is centred on the origin: x ∈ [−40, 40] nm here)
                set_region(mesh, 1, Box(center = (-20e-9, 0, 0), size = (40e-9, 80e-9, 10e-9)))
                set_region(mesh, 2, Box(center = (20e-9, 0, 0), size = (40e-9, 80e-9, 10e-9)))
                K1 = 4.5e4
                K2 = 5.5e4
                sm3 = make_sm(sim3; design = (Ku = region_map(1 => K1, 2 => K2),))
                regions = Array(sim3.mesh.regions)
                @test count(==(1), regions) > 0 && count(==(2), regions) > 0

                run_forward!(sm3)
                G = gradient!(sm3).Ku
                @test G isa Vector{Float64} && length(G) == sim3.n_total
                g_ad = sum(G[regions .== 1])

                fd_region(δ) = begin
                    set_design!(sm3; Ku = region_map(1 => K1 + δ, 2 => K2))
                    run_forward!(sm3)
                end
                δ = 50.0
                Jp = fd_region(+δ)
                Jm = fd_region(-δ)
                set_design!(sm3; Ku = region_map(1 => K1, 2 => K2))
                g_fd = (Jp - Jm) / (2δ)
                err = abs(g_ad - g_fd) / max(abs(g_fd), abs(g_ad), 1e-30)
                println(rpad("  ✓ region 1 aggregate", 34),
                        "  adjoint = ", @sprintf("%+.6e", g_ad),
                        "   fd = ", @sprintf("%+.6e", g_fd),
                        "   rel err = ", @sprintf("%.2e", err))
                @test err < 1e-5
            end

            # ---- vacuum strip: Ms = 0 for i > 12 ----
            # masked loss, vacuum-tolerant m*_validation and P_t pass-through
            @testset "vacuum strip" begin
                ms_fun(i, j, k, dx, dy, dz) = i <= 12 ? Ms0 : 0.0
                sim4 = base_sim(; ms = ms_fun)
                sm4 = make_sm(sim4; design = (Ku = Ku0,))
                J = run_forward!(sm4)
                # vacuum cells contribute nothing: |target| = |m*| = 0 there
                vac = 13:16
                m = reshape(sm4.m_star, 3, :)
                @test maximum(abs.(m[:, vac])) == 0.0
                kv = Ref(Ku0)
                fd_check("vacuum strip, scalar Ku", sm4, G -> G.Ku,
                         s -> (kv[] += s; set_design!(sm4; Ku = kv[])), 50.0)
            end
        end

        println("test_statematch_fd.jl done.")
    end
end
