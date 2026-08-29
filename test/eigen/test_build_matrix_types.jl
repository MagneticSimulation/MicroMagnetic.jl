using MicroMagnetic
using Test
using SparseArrays
using LinearAlgebra

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
    sim = Sim(mesh; driver="LLG", integrator="Heun")
    set_Ms(sim, 8.0e5)
    init_m0(sim, m0_dir)
    if add_exchange
        MicroMagnetic.add_exch(sim, 1.3e-11)
    end
    if add_dmi
        # D along z so that Dz≠0 — the 1x1x4 mesh only has ±z bonds, so
        # DMI's coupling to z-bonds is via Dz.  With D=(Dx,0,0) the kernel
        # would contribute 0 everywhere (Dz=0 → 0 * FlatTerm = 0::F64),
        # which gives a false negative on every DMI symbolic test.
        MicroMagnetic.add_dmi(sim, (0.0, 0.0, 1e-3))
    end
    # Inject one Epsilon on seed_site, seed_axis (1-indexed)
    if seed_site >= 1 && seed_axis in (1, 2, 3)
        N_sites = length(sim.spin) ÷ 3
        1 <= seed_site <= N_sites || error("seed_site out of range")
        idx = 3 * (seed_site - 1) + seed_axis
        # Epsilon has no type parameter — use the (id, value) constructor.
        ε_j = MicroMagnetic.Epsilon(idx, 1.0)
        # cpart() extracts the Float64 constant part without triggering the
        # loud convert stub (which now refuses FlatTerm/Epsilon → Float64).
        old = MicroMagnetic.cpart(sim.spin[idx])
        sim.spin[idx] = old + ε_j   # FlatTerm(old, (idx=>1.0)) equiv
    end
    return sim
end

function _classify_scalar(x)
    Tname = nameof(typeof(x))
    if Tname === :FlatTerm
        # A FlatTerm with empty coefs is effectively a plain Float64
        # constant (e.g. zero(FlatTerm) at a non-perturbed site, or
        # 0 * FlatTerm propagating through a kernel).  Only count it as
        # symbolic if it actually carries ε coefficients.
        return isempty(x.coefs) ? :F64 : :FlatTerm_nε_gt0
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
    MicroMagnetic.effective_field_local(sim, sim.spin)
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
        @test r.n_sym == 0          # all F64 (no ε anywhere)
    end
    let sim = _setup_sim_with_seed(seed_site=0, add_dmi=true)
        r = _probe(sim)
        @test r.n_sym == 0          # no seed → still all F64
    end

    # ── A2: seed on site 2, axis x, Exch-only ─────────────────────────
    # Axis x couples to Exchange cross terms: Exchange.field *must* carry
    # the FlatTerm(ε_{site2, mx=axis1}) contribution on at least 1 entry.
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=1, add_dmi=false)
        r = _probe(sim)
        @test r.n_sym >= 3
        ex_entry = only(p for p in r.per_interaction if contains(first(p), "Exchange"))
        @test last(ex_entry).n_sym >= 1
    end

    # ── B2: same seed, but WITH DMI added, seed on axis 1 (mx) ────────
    # For a z-only mesh with Dz≠0, BulkDMI's cross_y on a z-bond is
    #   cross_y(0,0,±1, mk_x, mk_y, mk_z) = ±1*mk_x - 0*mk_z = ±mk_x.
    # So H_DMI_y at site i DOES depend on m_x at the neighbour.  A pure
    # δm_x perturbation at site 2 therefore propagates ε into
    # BulkDMI.field at sites 1 and 3 (each on its y-component).
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=1, add_dmi=true)
        r = _probe(sim)
        dmi_entry = only(p for p in r.per_interaction if first(p) == "BulkDMI")
        @test last(dmi_entry).n_sym >= 1
    end

    # ── B3: seed on site 2, axis 3 (mz) + DMI (z-only mesh) ───────────
    # On a 1x1x4 mesh the only bonds are ±z.  For z-bonds the BulkDMI
    # cross products are:
    #   cross_x = ay*m_z - az*m_y = 0*m_z - 1*m_y = -m_y
    #   cross_y = az*m_x - ax*m_z = 1*m_x - 0*m_z =  m_x
    #   cross_z = 0
    # So H_DMI only couples to (m_x, m_y) of neighbours — a pure δm_z
    # perturbation does NOT enter H_DMI on this mesh.  This is a true
    # negative: it asserts that the symbolic kernel does NOT spuriously
    # invent ε where the physics says there is none.
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=3, add_dmi=true)
        r = _probe(sim)
        dmi_entry = only(p for p in r.per_interaction if first(p) == "BulkDMI")
        # DMI.field symbolic count must be zero (no coupling to mz on z-bonds)
        @test last(dmi_entry).n_sym == 0
        # Exchange DOES couple to mz: fz at sites 1, 2, 3 each carry ε.
        @test r.n_sym >= 3
    end

    # ── B4: seed on site 2, axis 2 (my) + DMI (z-only mesh) ───────────
    # my DOES couple into H_DMI for z-bonds (cross_x = -m_y), so this is
    # the positive counterpart to B3: BulkDMI.field must now carry
    # genuine ε contributions, and sim.field must have strictly more
    # symbolic entries than the Exch-only A2 baseline.
    let sim = _setup_sim_with_seed(seed_site=2, seed_axis=2, add_dmi=true)
        r = _probe(sim)
        dmi_entry = only(p for p in r.per_interaction if first(p) == "BulkDMI")
        @test last(dmi_entry).n_sym >= 1
        # Exchange contributes 3 (fy at sites 1,2,3) + DMI contributes ≥1
        @test r.n_sym >= 4
    end

    # ── Physics spot-check: m0 == ẑ opens DMI torque in B matrix ──────
    # On a z-only-bond mesh (1×1×nz), BulkDMI on z-bonds gives
    #   H_DMI_x = -m_y(neighbour), H_DMI_y = m_x(neighbour), H_DMI_z = 0.
    # The resulting torque m × H_DMI lies in the (x,y)-plane.  For this
    # torque to survive the tangent-frame projection it must be
    # PERPENDICULAR to m0 — so m0 must have NO in-plane component that
    # aligns with the torque.  Choosing m0 = ẑ makes the tangent plane
    # the (x,y)-plane, so the full DMI torque is captured.
    # (m0 = x̂ is the false-positive trap: the DMI torque is then along
    # x̂ ∥ m0 and is projected out, giving B1 == B0 even though DMI is
    # present — the kernel is correct, the physics just has no projection.)
    let
        function build_for(; use_dmi::Bool)
            MicroMagnetic.set_precision(AbstractFloat)
            mesh = FDMesh(nx=1, ny=1, nz=2, dx=5e-9, dy=5e-9, dz=5e-9)
            sim = Sim(mesh; driver="LLG", integrator="Heun")
            set_Ms(sim, 8.0e5)
            init_m0(sim, (0.0, 0.0, 1.0))   # m0 ∥ ẑ → DMI torque ⊥ m0
            MicroMagnetic.add_exch(sim, 1.3e-11)
            use_dmi && MicroMagnetic.add_dmi(sim, (0.0, 0.0, 5e-3))
            build_matrix(sim)
        end
        B0 = build_for(use_dmi=false)
        B1 = build_for(use_dmi=true)
        rel = opnorm(B1 - B0) / max(opnorm(B0), 1e-300)
        @test rel > 1e-3
    end

end
