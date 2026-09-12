# Tests for the adaptive mesh refinement module (src/amr), following
# García-Cervera & Roma, IEEE Trans. Magn. 42(6), 1648 (2006).
#
# The GPSM integrator of the package is CPU-only (sparse CHOLMOD), so these
# tests run on CPU with Float64 only.

using Test
using MicroMagnetic
import KernelAbstractions
using Logging

set_backend("cpu")
set_precision(Float64)

const _AMR_TMP = joinpath(tempdir(), "amr_tests")
isdir(_AMR_TMP) || mkpath(_AMR_TMP)
cd(_AMR_TMP)

"""
Berger-Rigoutsos clustering invariants: disjoint rectangles covering every
flagged cell.
"""
function test_br_cluster()
    nx, ny = 32, 16
    flags = zeros(Bool, nx * ny)
    # diagonal band of flags
    for i in 1:nx
        j = round(Int, 1 + (i - 1) * (ny - 1) / (nx - 1))
        for dj in -1:1
            jj = clamp(j + dj, 1, ny)
            flags[MicroMagnetic._cell_index(i, jj, 1, nx, ny)] = true
        end
    end
    nflags = count(flags)
    rects = MicroMagnetic._br_cluster(flags, nx, ny, 1,
                                      MicroMagnetic.AMRRect(1:nx, 1:ny, 1:1), 0.7)
    @test length(rects) > 0
    cover = zeros(Bool, nx * ny)
    nflag_in_rects = 0
    for r in rects, jj in r.jr, ii in r.ir
        I = MicroMagnetic._cell_index(ii, jj, 1, nx, ny)
        @test !cover[I]  # rectangles must be disjoint
        cover[I] = true
        flags[I] && (nflag_in_rects += 1)
    end
    @test nflag_in_rects == nflags  # all flagged cells must be covered
    # refinement ratio: patches built from doubled flag rects stay even-aligned
    finer = MicroMagnetic._rects_to_finer(rects, (true, true, false))
    @test all(first(r.ir) % 2 == 1 && last(r.ir) % 2 == 0 for r in finer)
end

"""
The cell-centered coarse->fine interpolation reproduces linear fields
exactly and the fine->coarse restriction averages children.
"""
function test_interp_restrict()
    nx0, ny0 = 8, 4
    # linear magnetization field in the cell-center coordinates
    mfun(i, j, dx, dy) = (0.1 + 0.3 * (i - 0.5) * dx, -0.2 + 0.1 * (j - 0.5) * dy, 0.9)
    src = zeros(3 * nx0 * ny0)
    for j in 1:ny0, i in 1:nx0
        I = MicroMagnetic._cell_index(i, j, 1, nx0, ny0)
        mx, my, mz = mfun(i, j, 1.0, 1.0)
        src[3I-2] = mx; src[3I-1] = my; src[3I] = mz
    end
    dst = zeros(3 * nx0 * 2 * ny0 * 2)
    kernel! = MicroMagnetic.amr_interp_box_kernel!(KernelAbstractions.CPU(),
                                                   MicroMagnetic.groupsize[])
    # with normalize=false the operator reproduces a linear field exactly
    kernel!(dst, src, true, true, false, 2nx0, 2ny0, 1, nx0, ny0, 1, false;
            ndrange=(2nx0, 2ny0, 1))
    err = 0.0
    # skip the outermost fine cells, whose stencil is constant-extended
    for j in 2:(2ny0 - 1), i in 2:(2nx0 - 1)
        I = MicroMagnetic._cell_index(i, j, 1, 2nx0, 2ny0)
        mx, my, mz = mfun(i, j, 0.5, 0.5)
        err = max(err, abs(dst[3I-2] - mx), abs(dst[3I-1] - my), abs(dst[3I] - mz))
    end
    @test err < 1e-13
    # with normalize=true the interpolated magnetization is a unit field
    kernel!(dst, src, true, true, false, 2nx0, 2ny0, 1, nx0, ny0, 1, true;
            ndrange=(2nx0, 2ny0, 1))
    nerr = 0.0
    for I in 1:(2nx0 * 2ny0)
        nrm = sqrt(dst[3I-2]^2 + dst[3I-1]^2 + dst[3I]^2)
        nerr = max(nerr, abs(nrm - 1.0))
    end
    @test nerr < 1e-12

    # restriction: a constant unit-magnetization field is reproduced
    fine = zeros(3 * nx0 * 2 * ny0 * 2)
    fine[1:3:end] .= 1.0
    coarse = zeros(3 * nx0 * ny0)
    covered = ones(Bool, nx0 * ny0)
    kernel! = MicroMagnetic.amr_restrict_kernel!(KernelAbstractions.CPU(),
                                                 MicroMagnetic.groupsize[])
    kernel!(coarse, fine, true, true, false, nx0, ny0, covered, true;
            ndrange=(nx0, ny0, 1))
    @test maximum(abs, coarse[1:3:end] .- 1.0) < 1e-14
    @test maximum(abs, coarse[2:3:end]) < 1e-14
    @test maximum(abs, coarse[3:3:end]) < 1e-14
end

"""
AMRSim with a single level reproduces the package's uniform GPSM integrator
step by step (same formulas, same dt on the whole grid).
"""
function test_gpsm_equivalence()
    mesh = FDMesh(dx=5e-9, dy=5e-9, nx=16, ny=8, nz=1)

    sim = Sim(mesh; name="_amr_eq", integrator="GPSM", save_data=false)
    set_Ms(sim, 8e5)
    add_exch(sim, 1.3e-11)
    add_anis(sim, 1e5; axis=(0, 0, 1))
    set_alpha(sim, 0.1)
    init_m0(sim, (1, 0.1, 0))
    sim.driver.integrator.step = 1e-13

    amr = AMRSim(mesh; levels=1, Ms=8e5, A=1.3e-11, Ku=1e5, alpha=0.1,
                 demag=false, name="_amr_eq1", save_data=false)
    init_m0(amr, (1, 0.1, 0))

    dt = 1e-13
    for step in 1:10
        MicroMagnetic.advance_step(sim)
        MicroMagnetic.amr_step!(amr, dt)
    end
    a = Array(sim.spin)
    b = Array(amr.base.sim.spin)
    @test maximum(abs.(a .- b)) < 1e-10
end

"""
A refined composite grid (exchange + anisotropy relaxation) converges to the
same relaxed state as the equivalent uniform finest grid.
"""
function test_amr_vs_uniform_relax()
    ms = 8e5
    Aex = 1.3e-11
    h0 = 2e5
    # 180 deg domain wall across the middle of the sample: the divergence
    # flags localize on the wall, a +x Zeeman field drives it out, and the
    # unique ground state m = +x is grid independent
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        mx = 40e-9 <= x <= 60e-9 ? -1.0 : 1.0
        return (mx, 0.1, 0.0)
    end

    # composite: 32x8 base grid + 1 refinement level (= 64x16 finest)
    mesh = FDMesh(dx=5e-9, dy=5e-9, nx=32, ny=8, nz=1)
    amr = AMRSim(mesh; levels=2, Ms=ms, A=Aex, H0=(h0, 0.0, 0.0), alpha=0.1,
                 demag=false, remesh_interval=50,
                 name="_amr_relax", save_data=false)
    init_m0(amr, m0fun)
    @test sum(length.(amr.patches)) > 0  # wall flags must generate patches
    relax(amr; dt=1e-12, maxsteps=4000, stopping_dmdt=1e-3, verbose=false)

    # uniform reference on the finest resolution
    mesh_u = FDMesh(dx=2.5e-9, dy=2.5e-9, nx=64, ny=16, nz=1)
    sim = Sim(mesh_u; name="_amr_ref", integrator="GPSM", save_data=false)
    set_Ms(sim, ms)
    add_exch(sim, Aex)
    add_zeeman(sim, (h0, 0, 0))
    set_alpha(sim, 0.1)
    init_m0(sim, m0fun)
    sim.driver.integrator.step = 1e-12
    for i in 1:4000
        MicroMagnetic.advance_step(sim)
    end

    am = MicroMagnetic.average_m(amr)
    bm = MicroMagnetic.average_m(sim)
    @test maximum(abs.(collect(am) .- collect(bm))) < 0.01
    @test am[1] > 0.95  # both systems relax to the +x ground state
    @test all(isfinite, Array(amr.base.sim.spin))
end

"""
The layered demag boxes reproduce the uniform demag field: exact on the base
level and within interpolation accuracy inside refined patches.
"""
function test_amr_demag_field()
    ms = 8e5
    # Bloch wall, delta = 5 nm: the user's baseline methodology -- restrict the
    # finest-grid field to every level's cells and compare (Error-A floor =
    # the uniform-same-resolution error).
    m0fun(i, j, k, dx, dy, dz) = begin
        x = (i - 0.5) * dx
        s = (x - 40e-9) / 5e-9
        (tanh(s), 0.0, sech(s))
    end
    function fine_ref(nfx)
        mesh = FDMesh(dx=1.25e-9, dy=1.25e-9, nx=nfx, ny=nfx ÷ 4, nz=1)
        sim = Sim(mesh; name="_amr_demag_ref$nfx", save_data=false)
        set_Ms(sim, ms)
        add_demag(sim)
        init_m0(sim, m0fun)
        demag = sim.interactions[findfirst(x -> isa(x, MicroMagnetic.Demag),
                                           sim.interactions)]
        MicroMagnetic.effective_field(demag, sim, sim.spin, 0.0)
        return Array(demag.field)
    end
    function restrict(href, nfx, ncx, ratio)
        nc = (ncx, ncx ÷ 4)
        out = zeros(3 * nc[1] * nc[2])
        for jj in 1:nc[2], ii in 1:nc[1]
            Ic = MicroMagnetic._cell_index(ii, jj, 1, nc[1], nc[2])
            acc = zeros(3)
            for dj in 0:ratio-1, di in 0:ratio-1
                If = MicroMagnetic._cell_index((ii - 1) * ratio + 1 + di,
                                               (jj - 1) * ratio + 1 + dj, 1, nfx, nfx ÷ 4)
                for d in 1:3
                    acc[d] += href[3If - 3 + d]
                end
            end
            for d in 1:3
                out[3Ic - 3 + d] = acc[d] / ratio^2
            end
        end
        return out
    end
    function level_err(amr, l, href, nfx)
        ncx = amr.dims[l][1]
        tgt = restrict(href, nfx, ncx, nfx ÷ ncx)
        inpatch = l == 1 ? zeros(Bool, length(amr.covered[1])) : amr.auth[l]
        cov = amr.covered[l]
        P = Array(amr.Phi[l])
        e = 0.0
        for i in 1:length(inpatch)
            (inpatch[i] && !cov[i]) || continue   # authoritative patch cells
            for d in 1:3
                e = max(e, abs(P[3i - 3 + d] - tgt[3i - 3 + d]))
            end
        end
        return e
    end
    function base_context_err(amr, href, nfx)
        tgt = restrict(href, nfx, amr.dims[1][1], nfx ÷ amr.dims[1][1])
        cov = amr.covered[1]
        P = Array(amr.Phi[1])
        e = 0.0
        for i in 1:length(cov)
            cov[i] && continue   # authoritative base cells
            for d in 1:3
                e = max(e, abs(P[3i - 3 + d] - tgt[3i - 3 + d]))
            end
        end
        return e
    end

    # config A: full coverage (patches refine the whole domain)
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=8, ny=2, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=ms, A=1.3e-11, demag=true,
                 remesh_interval=1000, name="_amr_demagA", save_data=false)
    init_m0(amr, m0fun)
    MicroMagnetic.sync_composite!(amr)
    MicroMagnetic.update_demag_boxes!(amr)
    href64 = fine_ref(64)
    hmax = maximum(abs, href64)
    e3 = level_err(amr, 3, href64, 64) / hmax
    e1 = base_context_err(amr, href64, 64) / hmax
    @test e3 < 0.02   # fine patch cells: at the uniform-same-resolution floor
    @test e1 < 0.02   # base cells: the composite beats the uniform coarse grid

    # config B: partial coverage (local refinement around the wall)
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=16, ny=4, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=ms, A=1.3e-11, demag=true,
                 remesh_interval=1000, box_buffer=8, refine_boundary=false,
                 name="_amr_demagB", save_data=false)
    init_m0(amr, m0fun)
    MicroMagnetic.sync_composite!(amr)
    MicroMagnetic.update_demag_boxes!(amr)
    href128 = fine_ref(128)
    hmax128 = maximum(abs, href128)
    e3b = level_err(amr, 3, href128, 128) / hmax128
    e1b = base_context_err(amr, href128, 128) / hmax128
    @test length(amr.patches[1]) > 0 && length(amr.patches[2]) > 0
    @test e3b < 0.02
    @test e1b < 0.01
    @test e1b < 0.04   # clearly better than the uniform 10 nm grid (0.040)
end

"""
Paper-style smoke run (Sections II/IV): Permalloy-like rectangle with strong
out-of-plane anisotropy, demag on, several refinement levels, picosecond
steps. Checks stability, energy dissipation and bounded refinement coverage.
"""
function test_diamond_smoke()
    mesh = FDMesh(dx=5e-9, dy=5e-9, nx=48, ny=12, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=8e5, A=1.3e-11, Ku=1e4, alpha=0.02,
                 H0=(2e5, 0.0, 0.0), demag=true, dt=5e-13, remesh_interval=50,
                 name="_amr_diamond", save_data=false)
    # 180 degree wall: demag charges on the wall + zeeman drive -> energy decreases
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        mx = 100e-9 <= x <= 140e-9 ? -1.0 : 1.0
        return (mx, 0.1, 0.0)
    end
    init_m0(amr, m0fun)
    MicroMagnetic.amr_step!(amr, 5e-13)   # builds the demag boxes for a consistent e0
    e0 = MicroMagnetic.amr_energies(amr).total
    relax(amr; dt=5e-13, maxsteps=300, stopping_dmdt=-1.0, verbose=false)
    e1 = MicroMagnetic.amr_energies(amr).total
    spin = Array(amr.base.sim.spin)
    @test all(isfinite, spin)
    @test e1 < e0  # the relaxed energy must decrease
    # refinement coverage stays bounded
    cov = count(amr.covered[2]) / length(amr.covered[2])
    @test cov <= 0.9
    save_ovf(amr, joinpath(_AMR_TMP, "_amr_diamond"))
    @test isfile(joinpath(_AMR_TMP, "_amr_diamond.ovf"))
end


function test_remesh_nesting()
    # repeated remeshing must keep the patch hierarchy properly nested and the
    # composite coverage telescoping (every physical cell claimed exactly once)
    ms = 8e5
    m0fun(i, j, k, dx, dy, dz) = begin
        x = (i - 0.5) * dx
        s = (x - 40e-9) / 5e-9
        (tanh(s), 0.0, sech(s))
    end
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=32, ny=8, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=ms, A=1.3e-11, Ku=5e4, demag=true,
                 alpha=0.02, remesh_interval=2, name="_amr_nest", save_data=false)
    init_m0(amr, m0fun)
    Vtrue = 32 * 10e-9 * 8 * 10e-9 * 1e-9
    for step in 1:6
        MicroMagnetic.amr_step!(amr, 5e-13)
        # proper nesting: every level-(l+1) patch lies inside a level-l patch
        for l in 3:amr.levels
            for p in amr.patches[l - 1]
                pir = ceil(Int, p.ir[1] / 2):ceil(Int, p.ir[end] / 2)
                pjr = ceil(Int, p.jr[1] / 2):ceil(Int, p.jr[end] / 2)
                pkr = ceil(Int, p.kr[1] / 2):ceil(Int, p.kr[end] / 2)
                ok = any(q.ir[1] <= pir[1] && pir[end] <= q.ir[end] &&
                         q.jr[1] <= pjr[1] && pjr[end] <= q.jr[end] &&
                         q.kr[1] <= pkr[1] && pkr[end] <= q.kr[end] for q in amr.patches[l - 2])
                @test ok
            end
        end
        # telescoping coverage: base non-covered + patch interiors non-covered
        # (each at its own level volume) == physical sample volume
        vol = 0.0
        for l in 1:amr.levels
            Vl = MicroMagnetic._level_volume(amr, l)
            if l == 1
                cnt = sum(!amr.covered[1][i] for i in 1:length(amr.covered[1]))
                vol += cnt * Vl
            else
                for p in amr.patches[l - 1]
                    cov = amr.covered[l]
                    for c in 1:length(p.kr), b in 1:length(p.jr), a in 1:length(p.ir)
                        Ig = MicroMagnetic._cell_index(first(p.ir) - 1 + a,
                                                       first(p.jr) - 1 + b,
                                                       first(p.kr) - 1 + c,
                                                       amr.dims[l][1], amr.dims[l][2])
                        cov[Ig] && continue
                        vol += Vl
                    end
                end
            end
        end
        @test abs(vol - Vtrue) / Vtrue < 1e-12
    end
end

"""
The saver's five energy columns share one `amr_energies` evaluation per saved
row: the cache returns the same object until the next step invalidates it, its
values match a fresh evaluation field by field, and the saved E_total column
matches a fresh evaluation of the same state.
"""
function test_saver_energy_cache()
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        mx = 100e-9 <= x <= 140e-9 ? -1.0 : 1.0
        return (mx, 0.1, 0.0)
    end
    mesh = FDMesh(dx=5e-9, dy=5e-9, nx=48, ny=12, nz=1)
    amr = AMRSim(mesh; levels=2, Ms=8e5, A=1.3e-11, Ku=1e4, alpha=0.02,
                 H0=(2e5, 0.0, 0.0), demag=true, dt=5e-13, remesh_interval=50,
                 name="_amr_ecache", save_data=true)
    init_m0(amr, m0fun)
    @test amr.ecache === nothing          # init_m0 invalidates the cache
    MicroMagnetic.amr_step!(amr, 5e-13)
    cached = MicroMagnetic._amr_energies_cached(amr)
    @test cached === MicroMagnetic._amr_energies_cached(amr)   # cache hit
    fresh = MicroMagnetic.amr_energies(amr)
    for n in propertynames(fresh)
        @test getfield(cached, n) == getfield(fresh, n)        # bitwise equal
    end
    MicroMagnetic.amr_step!(amr, 5e-13)
    @test MicroMagnetic._amr_energies_cached(amr) !== cached   # refreshed next step
    # the saved E_total column matches a fresh evaluation of the same state
    rows = filter(r -> !startswith(r, "#"),
                  readlines(joinpath(_AMR_TMP, "_amr_ecache.txt")))
    @test length(rows) == 2
    etot_row = parse(Float64, split(rows[end])[3])
    etot = MicroMagnetic.amr_energies(amr).total
    @test abs(etot_row - etot) / abs(etot) < 1e-12
end

"""
init_m0 leaves the final initial composite grid in place: no extra remesh is
triggered at nsteps = 0 (for a function m0 the fine-structure re-flagging
remesh happens inside init_m0), and the cadence is `nsteps > 0 &&
nsteps % remesh_interval == 0`.
"""
function test_first_step_no_remesh()
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        s = (x - 40e-9) / 5e-9
        (tanh(s), 0.0, sech(s))
    end
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=32, ny=8, nz=1)
    # tuple m0: one remesh in init_m0, none at the first step
    amr = AMRSim(mesh; levels=2, Ms=8e5, A=1.3e-11, demag=false,
                 remesh_interval=2, name="_amr_firststep", save_data=false)
    init_m0(amr, (1.0, 0.1, 0.0))
    @test amr.n_remesh == 1
    MicroMagnetic.amr_step!(amr, 5e-13)
    @test amr.n_remesh == 1
    # function m0: the fine-structure remesh is part of init_m0 (two in total)
    amr2 = AMRSim(mesh; levels=2, Ms=8e5, A=1.3e-11, demag=false,
                  remesh_interval=2, name="_amr_firststep2", save_data=false)
    init_m0(amr2, m0fun)
    @test sum(length.(amr2.patches)) > 0
    @test amr2.n_remesh == 2
    MicroMagnetic.amr_step!(amr2, 5e-13)
    @test amr2.n_remesh == 2
    @test all(isfinite, Array(amr2.base.sim.spin))
    # remesh cadence with interval 2: at step start nsteps = 1 -> no, 2 -> yes
    MicroMagnetic.amr_step!(amr2, 5e-13)
    @test amr2.n_remesh == 2
    MicroMagnetic.amr_step!(amr2, 5e-13)
    @test amr2.n_remesh == 3
end

"""
Remeshing warns when the divergence-threshold doubling loop is capped and the
finest-level coverage still exceeds `max_coverage` (deterministically triggered
with max_coverage = 0: the boundary flags do not depend on the threshold).
"""
function test_remesh_coverage_warn()
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=16, ny=8, nz=1)
    @test_logs (:warn, r"max_coverage") match_mode = :any begin
        amr = AMRSim(mesh; levels=2, Ms=8e5, A=1.3e-11, demag=false,
                     max_coverage=0.0, remesh_interval=1000,
                     name="_amr_covwarn", save_data=false)
        init_m0(amr, (1.0, 0.1, 0.0))
    end
end

"""
tensor2box two-box Newell kernel checks: equal sizes reproduce the package's
27-point stencil (up to summation order), splitting a source box into
sub-boxes is additive (which pins the /V_target normalization: a /V_source
convention would break the identity by the splitting factor), and the far
field approaches the analytic dipole in sign, direction and magnitude.
"""
function test_t2b_kernel()
    t2b = MicroMagnetic.demag_tensor_2box
    t2bmat(x, y, z, s, t) = begin
        xx = t2b(1, x, y, z, s..., t...)
        yy = t2b(2, x, y, z, s..., t...)
        zz = t2b(3, x, y, z, s..., t...)
        xy = t2b(4, x, y, z, s..., t...)
        xz = t2b(5, x, y, z, s..., t...)
        yz = t2b(6, x, y, z, s..., t...)
        return [xx xy xz; xy yy yz; xz yz zz]
    end
    # (a) equal sizes vs the package stencil: the 64-corner accumulation is
    #     not term-by-term the 27-point formula, so equality holds to ~1e-15
    #     (last-ulp summation order), not bitwise; the atol floor covers
    #     components that are exactly zero in one formula and ~1e-17 in the other
    for (x, y, z) in ((0.0, 0.0, 0.0), (7.5, -3.0, 0.0), (12.0, 12.0, 5.0),
                      (-20.0, 4.0, 9.0))
        h = 5.0
        @test isapprox(t2b(1, x, y, z, h, h, h, h, h, h),
                       MicroMagnetic.demag_tensor_xx(x, y, z, h, h, h);
                       rtol=1e-13, atol=1e-14)
        @test isapprox(t2b(2, x, y, z, h, h, h, h, h, h),
                       MicroMagnetic.demag_tensor_yy(x, y, z, h, h, h);
                       rtol=1e-13, atol=1e-14)
        @test isapprox(t2b(3, x, y, z, h, h, h, h, h, h),
                       MicroMagnetic.demag_tensor_zz(x, y, z, h, h, h);
                       rtol=1e-13, atol=1e-14)
        @test isapprox(t2b(4, x, y, z, h, h, h, h, h, h),
                       MicroMagnetic.demag_tensor_xy(x, y, z, h, h, h);
                       rtol=1e-13, atol=1e-14)
        @test isapprox(t2b(5, x, y, z, h, h, h, h, h, h),
                       MicroMagnetic.demag_tensor_xz(x, y, z, h, h, h);
                       rtol=1e-13, atol=1e-14)
        @test isapprox(t2b(6, x, y, z, h, h, h, h, h, h),
                       MicroMagnetic.demag_tensor_yz(x, y, z, h, h, h);
                       rtol=1e-13, atol=1e-14)
    end
    # (b) sub-box additivity: a 10x10x5 nm source equals the sum over 2x2
    #     5 nm subcells (contact, overlapping and far offsets)
    for tgt in ((5.0, 5.0, 5.0), (2.5, 2.5, 5.0)),
        (x, y, z) in ((7.5, 2.5, 0.0), (5.0, 5.0, 0.0), (30.0, -12.5, 0.0))
        for comp in 1:6
            big = t2b(comp, x, y, z, 10.0, 10.0, 5.0, tgt...)
            sub = sum(t2b(comp, x + ox, y + oy, z, 5.0, 5.0, 5.0, tgt...)
                      for ox in (-2.5, 2.5), oy in (-2.5, 2.5))
            @test isapprox(big, sub; rtol=1e-12, atol=1e-13)
        end
    end
    # (b2) target additivity: the parent-cell pairing equals the average of
    #      the equal-size children pairings -- the identity the down-gather's
    #      block average of the box FFT field relies on
    for (x, y) in ((7.5, 2.5), (5.0, 5.0), (30.0, -12.5), (-18.0, 11.5))
        for comp in 1:6
            big = t2b(comp, x, y, 0.0, 5.0, 5.0, 5.0, 10.0, 10.0, 5.0)
            sub = 0.25 * sum(t2b(comp, x + ox, y + oy, 0.0, 5.0, 5.0, 5.0,
                                 5.0, 5.0, 5.0)
                             for ox in (-2.5, 2.5), oy in (-2.5, 2.5))
            @test abs(big - sub) < 1e-13
        end
    end
    # (c) far field: H = -N.M (the down-gather pairing, Ms = 1) approaches
    #     the dipole field of the source cell; r/h ~ 50 keeps the corner-sum
    #     cancellation noise and the dipole-truncation error both < 1e-5
    h = 5.0
    V = h^3
    for (m, r) in (((1.0, 0.0, 0.0), (250.0, 100.0, 50.0)),
                   ((0.0, 1.0, 0.0), (-100.0, 225.0, 50.0)))
        H = -(t2bmat(r..., (h, h, h), (h, h, h)) * collect(m))
        r2 = r[1]^2 + r[2]^2 + r[3]^2
        mr = m[1] * r[1] + m[2] * r[2] + m[3] * r[3]
        ha = [(3 * mr * r[c] - m[c] * r2) / r2^2.5 * V / (4pi) for c in 1:3]
        @test maximum(abs.(H .- ha)) / maximum(abs.(ha)) < 1e-4
    end
    # 1/r^3 magnitude scaling along one direction (the ~3e-4 deviation from 8
    # is the dipole-truncation error of the analytic reference, not the kernel)
    r1 = (150.0, 75.0, 37.5)
    n1 = t2bmat(r1..., (h, h, h), (h, h, h))[1, 1]
    n2 = t2bmat(2 * r1[1], 2 * r1[2], 2 * r1[3], (h, h, h), (h, h, h))[1, 1]
    @test abs(n1 / n2 - 8) < 1e-2
    # reciprocity of the /V_target-normalized tensor: N_ST * V_T = N_TS^T * V_S
    nst = t2bmat(17.5, 7.5, 2.5, (10.0, 10.0, 5.0), (5.0, 5.0, 5.0))
    nts = t2bmat(17.5, 7.5, 2.5, (5.0, 5.0, 5.0), (10.0, 10.0, 5.0))
    @test maximum(abs.(nst .* 125.0 .- transpose(nts) .* 500.0)) < 1e-12
end

"""
Remesh pump canary: with the exact down-gather the energy jump of a remesh
step stays at the layout-change floor; an inconsistent (cubic-interpolated)
box gather produces order-of-magnitude larger remesh-step energy jumps in
long runs. Statistic: the largest per-remesh-step E_total increment over
1800 steps (48x12 base, 3 levels, remesh_interval=25; ~40 s). Threshold =
6x the post-fix maximum (1.32e-20 J, measured 2026-09-05); the cubic gather
measured 15.4e-20 J on the identical configuration (verified red via a
temporary git stash of the gather).
"""
function test_remesh_pump_canary()
    mesh = FDMesh(dx=7.8125e-9, dy=31.25e-9, dz=5e-9, nx=48, ny=12, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=8e5, A=1.3e-11, Ku=5e4, alpha=0.02,
                 demag=true, dt=5e-13, remesh_interval=25, name="_amr_pump",
                 save_data=false)
    init_m0(amr, (1, 0.2, 0))
    e_prev = MicroMagnetic.amr_energies(amr).total
    max_de = 0.0
    for s in 1:1800
        remeshes = amr.nsteps > 0 && amr.nsteps % amr.remesh_interval == 0
        MicroMagnetic.amr_step!(amr, 5e-13)
        e = MicroMagnetic.amr_energies(amr).total
        remeshes && (max_de = max(max_de, e - e_prev))
        e_prev = e
    end
    @test max_de < 8e-20
end

"""
Quiet assembly through the official `Sim` factory: construction, init_m0,
several steps and at least one remesh (rebuilding patches and demag boxes)
must (a) leave no saver files in the working directory, (b) not bump the
`_n_sims` counter, (c) not log any "has been created"/"has been added" @info.
"""
function test_quiet_assembly_side_effects()
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        s = (x - 40e-9) / 5e-9
        (tanh(s), 0.0, sech(s))
    end
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=32, ny=8, nz=1)
    n0 = MicroMagnetic._n_sims[]
    before = Set(readdir("."))
    logger = TestLogger(min_level=Logging.Info)
    amr = Logging.with_logger(logger) do
        amr = AMRSim(mesh; levels=2, Ms=8e5, A=1.3e-11, demag=true,
                     remesh_interval=2, name="_amr_quiet", save_data=false)
        init_m0(amr, m0fun)
        for _ in 1:5
            MicroMagnetic.amr_step!(amr, 5e-13)
        end
        amr
    end
    # (b) region sims are not user-visible sims
    @test MicroMagnetic._n_sims[] == n0
    # (c) no creation/addition @info leaks from the NullLogger windows
    msgs = [rec.message for rec in logger.logs]
    @test !any(contains(msg, "has been created") for msg in msgs)
    @test !any(contains(msg, "has been added") for msg in msgs)
    # (a) no saver files appear in the working directory
    new = setdiff(Set(readdir(".")), before)
    @test isempty(filter(f -> endswith(f, ".txt") || endswith(f, ".ovf"), collect(new)))
    # the region sims really come from the official factory ("None" driver)
    @test amr.base.sim.driver_name == "None"
end

"""
`Sim()`'s default behavior is untouched by the `quiet` kwarg: the creation
`@info` is still logged and the `_n_sims` counter still advances; with
`quiet=true` nothing is logged at info level and the counter stays.
"""
function test_sim_default_behavior_unchanged()
    n0 = MicroMagnetic._n_sims[]
    @test_logs (:info, "MicroSim (FD) has been created.") match_mode = :any begin
        Sim(FDMesh(nx=2, ny=2, nz=1))
    end
    @test MicroMagnetic._n_sims[] == n0 + 1
    n1 = MicroMagnetic._n_sims[]
    logger = TestLogger(min_level=Logging.Info)
    Logging.with_logger(logger) do
        Sim(FDMesh(nx=2, ny=2, nz=1); driver="None", name="_amr_quiet_sim",
            save_data=false, quiet=true)
    end
    @test MicroMagnetic._n_sims[] == n1
    @test isempty(logger.logs)
end

"""
P1 bookkeeping tool: `for_each_authoritative_cell` visits exactly the
authoritative composite cells -- per level, the `auth ∧ ¬covered` mask, counted
and index-summed independently below -- with `Ib` inside the region's own box,
and `composite_ncells` equals the tool count.
"""
function test_authoritative_cell_tool()
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        s = (x - 40e-9) / 5e-9
        (tanh(s), 0.0, sech(s))
    end
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=32, ny=8, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=8e5, A=1.3e-11, Ku=5e4, demag=true,
                 alpha=0.02, remesh_interval=1000, name="_amr_fec", save_data=false)
    init_m0(amr, m0fun)
    @test sum(length.(amr.patches)) > 0
    # independent per-level reference over the auth ∧ ¬covered masks
    nref = zeros(Int, amr.levels)
    sigref = zeros(Int, amr.levels)
    for l in 1:amr.levels
        for I in 1:length(amr.auth[l])
            if amr.auth[l][I] && !amr.covered[l][I]
                nref[l] += 1
                sigref[l] += I
            end
        end
    end
    ntool = zeros(Int, amr.levels)
    sigtool = zeros(Int, amr.levels)
    MicroMagnetic.for_each_authoritative_cell(amr) do r, Ib, Ig
        ntool[r.level] += 1
        sigtool[r.level] += Ig
        1 <= Ib <= r.sim.n_total || error("Ib outside the region box")
    end
    @test ntool == nref
    @test sigtool == sigref
    @test sum(ntool) == MicroMagnetic.composite_ncells(amr)
    # the single-region traversal covers the same cells as the filtered
    # whole-composite traversal (a region may legitimately hold zero
    # authoritative cells, e.g. a patch fully covered by finer patches)
    for r in MicroMagnetic.regions(amr)
        nsingle = Ref(0)
        MicroMagnetic.for_each_authoritative_cell(amr, r) do _r, _Ib, _Ig
            nsingle[] += 1
        end
        njoint = Ref(0)
        MicroMagnetic.for_each_authoritative_cell(amr) do rr, _Ib, _Ig
            rr === r && (njoint[] += 1)
        end
        @test nsingle[] == njoint[]
    end
end

"""
P3 kernels vs inline host references (verbatim copies of the pre-kernelization
triple loops): bitwise-equal divergence arrays and covered masks on small
configs with mixed refined dimensions and boundary cells, plus the end-to-end
`_divergence_arrays` and remesh-time masks on a real composite sim.
"""
function test_divergence_covered_kernels()
    # standalone divergence kernel vs inline reference (mixed refined dims)
    for (rx, ry, rz) in ((true, true, false), (true, false, false), (true, true, true))
        nx = rx ? 12 : 7
        ny = ry ? 9 : 5
        nz = rz ? 6 : 1
        C = zeros(3 * nx * ny * nz)
        for k in 1:nz, j in 1:ny, i in 1:nx
            I = MicroMagnetic._cell_index(i, j, k, nx, ny)
            C[3I-2] = 0.5 * sin(0.7i + 0.3j + 0.11k)
            C[3I-1] = 0.4 * cos(0.5i + 0.6j + 0.13k)
            C[3I] = 0.3 * sin(0.9i + 0.2j + 0.17k)
        end
        hx, hy, hz = 1.7, 0.9, 2.3
        div = zeros(nx * ny * nz)
        kernel! = MicroMagnetic.amr_divergence_kernel!(KernelAbstractions.CPU(),
                                                       MicroMagnetic.groupsize[])
        kernel!(div, C, rx, ry, rz, nx, ny, nz, hx, hy, hz; ndrange=(nx, ny, nz))
        ref = zeros(nx * ny * nz)
        for k in 1:nz, j in 1:ny, i in 1:nx
            I = MicroMagnetic._cell_index(i, j, k, nx, ny)
            d = 0.0
            if rx && nx > 1
                ip = i == nx ? i : i + 1
                im = i == 1 ? i : i - 1
                d += (C[3 * MicroMagnetic._cell_index(ip, j, k, nx, ny) - 2] -
                      C[3 * MicroMagnetic._cell_index(im, j, k, nx, ny) - 2]) /
                     ((ip - im) * hx)
            end
            if ry && ny > 1
                jp = j == ny ? j : j + 1
                jm = j == 1 ? j : j - 1
                d += (C[3 * MicroMagnetic._cell_index(i, jp, k, nx, ny) - 1] -
                      C[3 * MicroMagnetic._cell_index(i, jm, k, nx, ny) - 1]) /
                     ((jp - jm) * hy)
            end
            if rz && nz > 1
                kp = k == nz ? k : k + 1
                km = k == 1 ? k : k - 1
                d += (C[3 * MicroMagnetic._cell_index(i, j, kp, nx, ny)] -
                      C[3 * MicroMagnetic._cell_index(i, j, km, nx, ny)]) / ((kp - km) * hz)
            end
            ref[I] = abs(d)
        end
        @test div == ref
    end
    # standalone covered-mask kernel vs inline reference
    for (rx, ry, rz) in ((true, true, false), (true, true, true))
        ncx = rx ? 6 : 5
        ncy = ry ? 7 : 4
        ncz = rz ? 3 : 1
        nfx = ncx * (rx ? 2 : 1)
        nfy = ncy * (ry ? 2 : 1)
        nfz = rz ? 2 * ncz : ncz
        fine = zeros(Bool, nfx * nfy * nfz)
        for k in 1:nfz, j in 1:nfy, i in 1:nfx
            fine[MicroMagnetic._cell_index(i, j, k, nfx, nfy)] = (i + 2j + 3k) % 3 != 0
        end
        cov = zeros(Bool, ncx * ncy * ncz)
        kernel! = MicroMagnetic.amr_covered_mask_kernel!(KernelAbstractions.CPU(),
                                                         MicroMagnetic.groupsize[])
        kernel!(cov, fine, rx, ry, rz, ncx, ncy, nfx, nfy; ndrange=(ncx, ncy, ncz))
        ref = zeros(Bool, ncx * ncy * ncz)
        for k in 1:ncz, j in 1:ncy, i in 1:ncx
            ok = true
            for kk in (rz ? 2 * k - 1 : k, rz ? 2 * k : k),
                jj in (ry ? 2 * j - 1 : j, ry ? 2 * j : j),
                ii in (rx ? 2 * i - 1 : i, rx ? 2 * i : i)
                ok &= fine[MicroMagnetic._cell_index(ii, jj, kk, nfx, nfy)]
            end
            ref[MicroMagnetic._cell_index(i, j, k, ncx, ncy)] = ok
        end
        @test cov == ref
    end
    # end-to-end on a real composite sim: _divergence_arrays and the
    # remesh-time covered masks against inline references
    function m0fun(i, j, k, dx, dy, dz)
        x = (i - 0.5) * dx
        s = (x - 40e-9) / 5e-9
        (tanh(s), 0.0, sech(s))
    end
    mesh = FDMesh(dx=10e-9, dy=10e-9, nx=16, ny=4, nz=1)
    amr = AMRSim(mesh; levels=3, Ms=8e5, A=1.3e-11, Ku=5e4, demag=false,
                 alpha=0.02, remesh_interval=1000, name="_amr_dvk", save_data=false)
    init_m0(amr, m0fun)
    @test sum(length.(amr.patches)) > 0
    divs, divmax = MicroMagnetic._divergence_arrays(amr)
    refmax = 0.0
    for l in 1:(amr.levels - 1)
        nx, ny, nz = amr.dims[l]
        hx, hy, hz = MicroMagnetic._level_h(amr.base_mesh, amr.refined, l)
        rx, ry, rz = amr.refined
        C = Array(amr.C[l])
        ref = zeros(nx * ny * nz)
        for k in 1:nz, j in 1:ny, i in 1:nx
            I = MicroMagnetic._cell_index(i, j, k, nx, ny)
            d = 0.0
            if rx && nx > 1
                ip = i == nx ? i : i + 1
                im = i == 1 ? i : i - 1
                d += (C[3 * MicroMagnetic._cell_index(ip, j, k, nx, ny) - 2] -
                      C[3 * MicroMagnetic._cell_index(im, j, k, nx, ny) - 2]) /
                     ((ip - im) * hx)
            end
            if ry && ny > 1
                jp = j == ny ? j : j + 1
                jm = j == 1 ? j : j - 1
                d += (C[3 * MicroMagnetic._cell_index(i, jp, k, nx, ny) - 1] -
                      C[3 * MicroMagnetic._cell_index(i, jm, k, nx, ny) - 1]) /
                     ((jp - jm) * hy)
            end
            if rz && nz > 1
                kp = k == nz ? k : k + 1
                km = k == 1 ? k : k - 1
                d += (C[3 * MicroMagnetic._cell_index(i, j, kp, nx, ny)] -
                      C[3 * MicroMagnetic._cell_index(i, j, km, nx, ny)]) / ((kp - km) * hz)
            end
            ref[I] = abs(d)
        end
        @test divs[l] == ref
        refmax = max(refmax, maximum(ref))
    end
    @test divmax == refmax
    for l in 1:(amr.levels - 1)
        ncx, ncy, ncz = amr.dims[l]
        rx, ry, rz = amr.refined
        fine = amr.auth[l + 1]
        nfx, nfy = ncx * (rx ? 2 : 1), ncy * (ry ? 2 : 1)
        ref = zeros(Bool, ncx * ncy * ncz)
        for k in 1:ncz, j in 1:ncy, i in 1:ncx
            ok = true
            for kk in (rz ? 2 * k - 1 : k, rz ? 2 * k : k),
                jj in (ry ? 2 * j - 1 : j, ry ? 2 * j : j),
                ii in (rx ? 2 * i - 1 : i, rx ? 2 * i : i)
                ok &= fine[MicroMagnetic._cell_index(ii, jj, kk, nfx, nfy)]
            end
            ref[MicroMagnetic._cell_index(i, j, k, ncx, ncy)] = ok
        end
        @test amr.covered[l] == ref
    end
    @test !any(amr.covered[amr.levels])
end

@testset "AMR" begin
    test_br_cluster()
    test_interp_restrict()
    test_t2b_kernel()
    test_gpsm_equivalence()
    test_amr_vs_uniform_relax()
    test_amr_demag_field()
    test_diamond_smoke()
    test_remesh_nesting()
    test_saver_energy_cache()
    test_first_step_no_remesh()
    test_remesh_coverage_warn()
    test_authoritative_cell_tool()
    test_divergence_covered_kernels()
    test_remesh_pump_canary()
    test_quiet_assembly_side_effects()
    test_sim_default_behavior_unchanged()
end
