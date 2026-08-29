using MicroMagnetic
using Test

"""
Verify that the build_matrix 2N-ε symbolic pass correctly propagates
Epsilon / AddExpr *through* the KA CPU kernels.

Bugs fixed by the companion kernel change (fx = m[i]-m[i]+T(0)):
  * When a kernel declares fx = T(0) with T=Float64 while m[i] contains
    Epsilon / AddExpr values (which is what build_matrix's symbolic pass
    feeds into effective_field_local), the KernelAbstractions CPU backend
    silently truncates the widened accumulator back to Float64, so the
    ε-coefficients never reach H_eff.
  * The same silent truncation happens on every `energy[I] = 0` /
    `h[i..i+2] = 0` assignment in Ms==0 branches, converting the output
    buffers AbstractFloat element-by-element back towards Float64.

Both bugs are detected by probing interaction.field / sim.field for the
presence of ε-tagged values (FlatTerm with nε>0, or AddExpr) after a
single-site m perturbation on top of set_precision(AbstractFloat).
"""

const _no_seed = (0, 0, 0.0)     # sentinel: do not inject ε

function _setup_sim_with_seed(; add_exchange::Bool=true,
                                add_dmi::Bool=false,
                                seed_site::Int=2,
                                seed_axis::Int=1,   # 1=mx, 2=my, 3=mz
                                m0_dir=(0.0, 0.0, 1.0))
    set_precision(AbstractFloat)
    mesh = FDMesh(nx=1, ny=1, nz=4, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh; integrator="RUNGE_KUTTA4")
    set_Ms(sim, 8.0e5)
    init_m0(sim, m0_dir)
    if add_exchange
        add_exchange(sim, 1.3e-11)
    end
    if add_dmi
        add_dmi(sim, (1e-3, 0.0, 0.0))
    end
    # Inject one Epsilon on seed_site, seed_axis (1-indexed)
    if seed_site >= 1 && seed_axis in (1, 2, 3)
        N_sites = length(sim.spin) ÷ 3
        1 <= seed_site <= N_sites || error("seed_site out of range")
        idx = 3 * (seed_site - 1) + seed_axis
        ε_j = MicroMagnetic.Epsilon{idx}()
        old = convert(Float64, sim.spin[idx])
        sim.spin[idx] = old + ε_j   # FlatTerm(old, (idx=>1.0)) equiv
    end
    return sim
end

function _classify_scalar(x)
    Tname = nameof(typeof(x))
    if Tname === :FlatTerm
        return :FlatTerm_nε_gt0   # symbolic, one ε
    elseif Tname === :AddExpr
        return :AddExpr           # symbolic, many terms
    elseif x isa Float64
        return :F64               # no symbol attached (truncation bug!)
    else
        return :other
    end
end

function _count_symbolic_entries(arr)
    n_sym = 0
    n_f64 = 0
    for v in arr
        c = _classify_scalar(v)
        if c === :F64
            n_f64 += 1
        else
            n_sym += 1
        end
    end
    return (; n_sym, n_f64)
end

function _probe(sim)
    # Drive the exact same effective_field_local() code path build_matrix uses
    MicroMagnetic.effective_field_local(sim, sim.spin, 0.0)
    N_total = length(sim.field)
    s_stat = _count_symbolic_entries(sim.field)
    per_interaction = Pair{String, Any}[]
    for inter in sim.interactions
        T = typeof(inter)
        push!(per_interaction, string(nameof(T)) => _count_symbolic_entries(inter.field))
    end
    return (; N_total, s_stat..., per_interaction)
end

@testset "build_matrix — kernel symbolic-type propagation (Exchange / DMI)" begin

    # ── A1 / B1: no seed ──────────────────────────────────────────────
    let sim = _setup_sim_with_seed(seed_site=0, add_dmi=false)
        r = _probe(sim)
        @test r.s_stat.n_sym == 0          # all F64 (no ε anywhere)
    end
    let sim = _setup_sim_with_seed(seed_site=0, add_dmi=true)
        r = _probe(sim)
        @test r.s_stat.n_sym == 0          # no seed → still all F64
    end

    # ── A2: seed on site 2, axis x, Exch-only ─────────────────────────
    # Axis x couples to Exchange cross terms: Exchange.field *must* carry
    # the FlatTerm(ε_{site2, mx=axis1}) contribution on at least 1 entry.
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=1, add_dmi=false)
        r = _probe(sim)
        @test r.s_stat.n_sym >= 3
        ex_entry = only(p for p in r.per_interaction if first(p) === :Exchange)
        @test last(ex_entry).n_sym >= 1
    end

    # ── B2: same seed, but WITH DMI added, seed on axis 1 (mx) ────────
    # DMI is H ∝ cross(∇, m) and cross only depends on (my, mz) components
    # of the gradient of m; a pure δm_x perturbation at a single cell
    # does NOT enter H_DMI linearly (only indirectly through cross which
    # needs m_y or m_z variation).  So DMI.field MUST NOT suddenly show
    # symbols when seed is on mx.
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=1, add_dmi=true)
        r = _probe(sim)
        dmi_entry = only(p for p in r.per_interaction if first(p) === :BulkDMI)
        # Exact same number of symbolic entries as A2 (DMI contributes 0
        # symbols for axis-x seed).
        @test last(dmi_entry).n_sym == 0
    end

    # ── B3: seed on site 2, axis 3 (mz) + DMI ─────────────────────────
    # mz variation DOES feed into H_DMI cross terms, so BulkDMI.field
    # itself must now show FlatTerm/AddExpr entries.  This is the binary
    # proof that the kernel widening trick actually works, and this test
    # would FAIL (0 / 0 / same count) on the buggy "= T(0),T(0),T(0)"
    # version of bulkdmi_kernel!.
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=3, add_dmi=true)
        r = _probe(sim)
        dmi_entry = only(p for p in r.per_interaction if first(p) === :BulkDMI)
        # DMI.field symbolic count must be strictly positive
        @test last(dmi_entry).n_sym >= 1
        # And sim.field as a whole must have strictly MORE symbolic entries
        # than the Exch-only A2 case (which had ≥3).
        @test r.s_stat.n_sym >= 5
    end

    # ── Physics spot-check: m0 == x̂ opens DMI torque in B matrix ──────
    # On m0 = (1,0,0), H_DMI has z-component which gives nonzero
    # m0 × H_DMI = x̂ × ẑ  cross term → the B matrix with/without DMI
    # must differ by more than 1e-3 in norm.
    let
        function build_for(; use_dmi::Bool)
            set_precision(AbstractFloat)
            mesh = FDMesh(nx=1, ny=1, nz=2, dx=5e-9, dy=5e-9, dz=5e-9)
            sim = Sim(mesh; integrator="RUNGE_KUTTA4")
            set_Ms(sim, 8.0e5)
            init_m0(sim, (1.0, 0.0, 0.0))   # m0 along x, required so
            add_exchange(sim, 1.3e-11)       # m0×H_DMI != 0
            use_dmi && add_dmi(sim, (5e-3, 0.0, 0.0))
            dropzeros!(build_matrix(sim))
        end
        B0 = Matrix(build_for(use_dmi=false))
        B1 = Matrix(build_for(use_dmi=true))
        rel = opnorm(B1 - B0) / max(opnorm(B0), 1e-300)
        @test rel > 1e-3
    end

end
