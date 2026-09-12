# Demag v3 for the composite grid: per-level replicated-parent solves + local
# correction boxes. NO magnetization interpolation, NO covered-field
# subtraction anywhere.
#
# Telescoping source identity. Let C[ℓ] be the composite magnetization array
# at level ℓ (authoritative patch data; interpolated values in cells covered
# by finer patches). Define, on each level-ℓ patch footprint F_ℓ,
#     Δ_ℓ = C_patch − rep_{ℓ←ℓ−1}(C[ℓ−1])        (zero outside F_ℓ),
# the fine data minus the parent-level composite replicated to h_ℓ (exact:
# piecewise-constant replication preserves the physical source). The
# composite source then telescopes EXACTLY for arbitrary nesting:
#     S_composite = C₁ + Δ₂ + Δ₃ + ... + Δ_L,
# because over F_ℓ the parent context plus Δ_ℓ is C_patch. Each Δ_ℓ is
# compact, and has zero net dipole moment whenever the restriction is exact.
#
# Fields (by linearity; every solve on its native lattice):
#  * Phi[ℓ] = H(rep_{ℓ←ℓ−1}(C[ℓ−1]))   -- one full-level "shadow" FFT: the
#    parent-level composite source evaluated at h_ℓ resolution
#    + the Δ_ℓ box corrections of the level-ℓ patches (free-space FFT on
#    bbox+buffer boxes, added on the native box lattice)
#    + the Δ_{ℓ+1} box fields gathered down from finer boxes (exact two-box
#    Newell pairing within their buffers).
#    Far contributions of compact zero-dipole Δ's are dropped: that is the
#    O(h²)-smooth coarse-representation error (Error A), shared by every
#    composite method.
#  * Phi[1] = H(C₁) (the base FFT) + the level-2 Δ box fields gathered at
#    base cells inside the box buffers. No subtraction is needed anywhere:
#    the −rep(C₁) part of Δ cancels the base source over the patch footprint
#    exactly at the SOURCE level, so H(C₁) + H(Δ) is the composite field
#    wherever the box fields reach.

# ---------------------------------------------------------------- kernels

"""
Replicate the parent-level composite `C[ℓ−1]` onto the level-ℓ grid: each
level cell takes its containing parent cell's value (exact: piecewise-constant
replication preserves the physical source). `pdx/pdy/pdz` = the parent dims;
refined dims have ratio 2, unrefined dims are identical.
"""
@kernel function amr_replicate_kernel!(dst, @Const(parent_spin), rx::Bool, ry::Bool,
                                       rz::Bool, ndx::Int, ndy::Int,
                                       pdx::Int, pdy::Int, pdz::Int)
    i, j, k = @index(Global, NTuple)
    ip = rx ? clamp(((i - 1) >> 1) + 1, 1, pdx) : min(i, pdx)
    jp = ry ? clamp(((j - 1) >> 1) + 1, 1, pdy) : min(j, pdy)
    kp = rz ? clamp(((k - 1) >> 1) + 1, 1, pdz) : min(k, pdz)
    Ib = _cell_index(i, j, k, ndx, ndy)
    e = 3 * _cell_index(ip, jp, kp, pdx, pdy) - 2
    d = 3 * Ib - 2
    @inbounds begin
        dst[d] = parent_spin[e]
        dst[d + 1] = parent_spin[e + 1]
        dst[d + 2] = parent_spin[e + 2]
    end
end

"""
Scatter the Δ source of one patch into its box:
    Δ = patch.spin − R_ℓ(parent cells)
on patch cells (zero elsewhere; the caller zero-fills the box array). The
gather reads the level-ℓ replicated-base box, which is piecewise-constant at
the base scale, so the gather returns the parent's value exactly.

Index bookkeeping: `bx,by,bz` = the patch's first cell relative to the box
origin (box cells); `px,py,pz` = the patch's first cell in level-ℓ global
coords; `pgx,pgy,pgz` = the patch's ghost-layer counts (its arrays are
ghosted); `pnx,pny` = the patch's ghosted dims.
"""
@kernel function amr_delta_kernel!(box_src, @Const(patch_spin), @Const(rep_spin),
                                   bx::Int, by::Int, bz::Int,
                                   px::Int, py::Int, pz::Int,
                                   ngx::Int, ngy::Int,
                                   pgx::Int, pgy::Int, pgz::Int,
                                   pnx::Int, pny::Int,
                                   nbx::Int, nby::Int, nbz::Int)
    a, b, c = @index(Global, NTuple)
    Ib = _cell_index(a + bx, b + by, c + bz, ngx, ngy)
    gi = px - 1 + a
    gj = py - 1 + b
    gk = pz - 1 + c
    Ir = _cell_index(gi, gj, gk, nbx, nby)
    Ip = _cell_index(a + pgx, b + pgy, c + pgz, pnx, pny)
    d = 3 * Ib - 2
    e = 3 * Ir - 2
    dp = 3 * Ip - 2
    @inbounds begin
        box_src[d] = patch_spin[dp] - rep_spin[e]
        box_src[d + 1] = patch_spin[dp + 1] - rep_spin[e + 1]
        box_src[d + 2] = patch_spin[dp + 2] - rep_spin[e + 2]
    end
end

"""
Own-level box add: the box's Δ field at its native fine cells. The box FFT
solve IS the native-lattice field, so this is a cell-to-cell add (no
interpolation; sh == 1 in the old gather kernel).
"""
@kernel function amr_box_own_add_kernel!(dst, @Const(box_field),
                                         i0::Int, j0::Int, k0::Int,
                                         ngx::Int, ngy::Int, dnx::Int, dny::Int)
    i, j, k = @index(Global, NTuple)
    Ib = _cell_index(i0 - 1 + i, j0 - 1 + j, k0 - 1 + k, dnx, dny)
    Ic = _cell_index(i, j, k, ngx, ngy)
    d = 3 * Ib - 2
    e = 3 * Ic - 2
    @inbounds begin
        dst[d] += box_field[e]
        dst[d + 1] += box_field[e + 1]
        dst[d + 2] += box_field[e + 2]
    end
end

# ---------------------------------------------------------------- box manager

"""Expand a patch rectangle by `buffer` cells (clamped to the level domain)
on refined dimensions."""
function _expand_rect(ir, jr, kr, buffer, dims, refined)
    irx = refined[1] ? (max(1, ir[1] - buffer):min(dims[1], ir[end] + buffer)) : ir
    jry = refined[2] ? (max(1, jr[1] - buffer):min(dims[2], jr[end] + buffer)) : jr
    krz = refined[3] ? (max(1, kr[1] - buffer):min(dims[3], kr[end] + buffer)) : kr
    return AMRRect(irx, jry, krz)
end

"""Merge overlapping rectangles (greedy pass until fixpoint)."""
function _merge_rects(rects::Vector{AMRRect})
    out = copy(rects)
    changed = true
    while changed && length(out) > 1
        changed = false
        for i in 1:length(out), j in (i+1):length(out)
            overlap = out[i].ir[1] <= out[j].ir[end] && out[j].ir[1] <= out[i].ir[end] &&
                      out[i].jr[1] <= out[j].jr[end] && out[j].jr[1] <= out[i].jr[end] &&
                      out[i].kr[1] <= out[j].kr[end] && out[j].kr[1] <= out[i].kr[end]
            if overlap
                m = AMRRect(min(out[i].ir[1], out[j].ir[1]):max(out[i].ir[end], out[j].ir[end]),
                            min(out[i].jr[1], out[j].jr[1]):max(out[i].jr[end], out[j].jr[end]),
                            min(out[i].kr[1], out[j].kr[1]):max(out[i].kr[end], out[j].kr[end]))
                deleteat!(out, j)
                out[i] = m
                changed = true
                break
            end
        end
    end
    return out
end

"""(Re)build the per-level correction boxes from the current patches."""
function rebuild_boxes!(amr::AMRSim)
    empty!(amr.dboxes)
    for l in 2:amr.levels
        boxes = DemagBox[]
        if !isempty(amr.patches[l - 1])
            expanded = [_expand_rect(p.ir, p.jr, p.kr, amr.box_buffer, amr.dims[l],
                                     amr.refined) for p in amr.patches[l - 1]]
            merged = _merge_rects(expanded)
            for m in merged
                nxi, nyi, nzi = length(m.ir), length(m.jr), length(m.kr)
                hx, hy, hz = _level_h(amr.base_mesh, amr.refined, l)
                mesh = FDMesh(dx=hx, dy=hy, dz=hz, nx=nxi, ny=nyi, nz=nzi,
                              x0=amr.base_mesh.x0 + (m.ir[1] - 1) * hx,
                              y0=amr.base_mesh.y0 + (m.jr[1] - 1) * hy,
                              z0=amr.base_mesh.z0 + (m.kr[1] - 1) * hz)
                sim = _region_sim(mesh, amr.Ms, "$(amr.name)_box_L$l")
                _add_region_demag!(sim)
                members = [pi for (pi, p) in enumerate(amr.patches[l - 1])
                           if p.ir ⊆ m.ir && p.jr ⊆ m.jr && p.kr ⊆ m.kr]
                push!(boxes, DemagBox(l, m.ir, m.jr, m.kr, sim, members))
            end
        end
        push!(amr.dboxes, boxes)
    end
    return nothing
end

"""
Down-gather a box correction field onto the immediate parent level: every
parent cell takes the plain average of the box field over its 2×2×1 children
(refined dims; identity on unrefined dims). The box FFT field at a fine cell
IS the exact native-lattice volume-averaged pairing (the package's equal-size
Newell stencil), and volume averaging is additive over children, so this
block average is exactly the two-box pairing of the parent cell -- no
interpolation, no approximation (the exact-transfer replacement of the old
cubic Lagrange gather, whose O(1) error at Δ charge sheets pumped energy at
every remesh).
"""
@kernel function amr_box_to_parent_kernel!(dst, @Const(box_field),
                                           bx0::Int, by0::Int, bz0::Int,
                                           ngx::Int, ngy::Int,
                                           dnx::Int, dny::Int,
                                           ci0::Int, cj0::Int, ck0::Int,
                                           rx::Bool, ry::Bool, rz::Bool)
    i, j, k = @index(Global, NTuple)
    ip, jp, kp = ci0 - 1 + i, cj0 - 1 + j, ck0 - 1 + k
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    @inbounds for kk in (rz ? (2 * kp - 1:2 * kp) : kp),
                   jj in (ry ? (2 * jp - 1:2 * jp) : jp),
                   ii in (rx ? (2 * ip - 1:2 * ip) : ip)
        e = 3 * _cell_index(ii - bx0 + 1, jj - by0 + 1, kk - bz0 + 1, ngx, ngy) - 2
        s1 += box_field[e]
        s2 += box_field[e + 1]
        s3 += box_field[e + 2]
    end
    n = 1.0 / ((rx ? 2.0 : 1.0) * (ry ? 2.0 : 1.0) * (rz ? 2.0 : 1.0))
    Ic = _cell_index(ip, jp, kp, dnx, dny)
    d = 3 * Ic - 2
    @inbounds begin
        dst[d] += n * s1
        dst[d + 1] += n * s2
        dst[d + 2] += n * s3
    end
end

# ---------------------------------------------------------------- field update

"""
    update_demag_boxes!(amr)

Assemble `amr.Phi[1..levels]`, the demag field at every authoritative cell
(see the module header for the telescoping identity):

  * `Phi[1]` = the base FFT field of `C₁` plus the level-2 Δ box fields
    gathered at base cells inside the box buffers;
  * `Phi[ℓ]` = the full-level shadow solve of `rep(C[ℓ−1])` (the parent
    composite context at h_ℓ resolution) plus the Δ_ℓ box corrections (native
    lattice inside each box) plus the finer Δ box fields gathered down
    (exact two-box Newell pairing within their buffers).
"""
function update_demag_boxes!(amr::AMRSim)
    amr.demag || return nothing
    rx, ry, rz = amr.refined
    base = amr.base

    # (a) base solve: H(C₁). Phi[1] aliases the base Demag field buffer, so
    #     this must run before every box gather below.
    dem = base.sim.interactions[findfirst(x -> isa(x, Demag), base.sim.interactions)]
    effective_field(dem, base.sim, base.sim.spin, 0.0)

    # (b) per level: parent-context shadow solve on the full level box
    #     R_ℓ = rep_{ℓ←ℓ−1}(C[ℓ−1]); H(R_ℓ) = the level-(ℓ−1) composite source
    #     evaluated on the h_ℓ lattice (exact on native cells)
    for l in 2:amr.levels
        sh = amr.shadow[l - 1]
        ndx, ndy, ndz = amr.dims[l]
        pdx, pdy, pdz = amr.dims[l - 1]
        kernel! = amr_replicate_kernel!(get_backend(sh.spin), groupsize[])
        kernel!(sh.spin, amr.C[l - 1], rx, ry, rz, ndx, ndy, pdx, pdy, pdz;
                ndrange=(ndx, ndy, ndz))
        effective_field(sh, sh.spin, 0.0)
        copyto!(amr.Phi[l], sh.field)
    end

    # (c) per box: delta source (patch spin − replicated parent) -> box FFT
    for l in 2:amr.levels, box in amr.dboxes[l - 1]
        ngx, ngy, ngz = box.sim.mesh.nx, box.sim.mesh.ny, box.sim.mesh.nz
        fill!(box.sim.spin, 0.0)
        rep = amr.shadow[l - 1].spin
        for pi in box.patches
            p = amr.patches[l - 1][pi]
            nxi, nyi, nzi = length(p.ir), length(p.jr), length(p.kr)
            gx, gy, gz = _ghosts(amr)
            kernel! = amr_delta_kernel!(get_backend(box.sim.spin), groupsize[])
            kernel!(box.sim.spin, p.sim.spin, rep,
                    p.ir[1] - box.ir[1], p.jr[1] - box.jr[1], p.kr[1] - box.kr[1],
                    p.ir[1], p.jr[1], p.kr[1],
                    ngx, ngy, gx, gy, gz,
                    p.sim.mesh.nx, p.sim.mesh.ny,
                    amr.dims[l][1], amr.dims[l][2], amr.dims[l][3];
                    ndrange=(nxi, nyi, nzi))
        end
        bdem = box.sim.interactions[findfirst(x -> isa(x, Demag), box.sim.interactions)]
        effective_field(bdem, box.sim, box.sim.spin, 0.0)
    end

    # (d) own-level box adds: each box's Δ field at its own-level cells
    for l in 2:amr.levels, box in amr.dboxes[l - 1]
        ngx, ngy = box.sim.mesh.nx, box.sim.mesh.ny
        bfield = box.sim.interactions[findfirst(x -> isa(x, Demag),
                                                box.sim.interactions)].field
        own! = amr_box_own_add_kernel!(get_backend(amr.Phi[l]), groupsize[])
        own!(amr.Phi[l], bfield, box.ir[1], box.jr[1], box.kr[1],
             ngx, ngy, amr.dims[l][1], amr.dims[l][2];
             ndrange=(length(box.ir), length(box.jr), length(box.kr)))
    end

    # (e) down-gathers: each box's Δ field at its immediate-parent cells,
    #     exact two-box pairing via the block average of the native-lattice
    #     box field (parent cells under the patch are covered/
    #     non-authoritative, so adding there is harmless)
    for l in 2:amr.levels, box in amr.dboxes[l - 1]
        ngx, ngy = box.sim.mesh.nx, box.sim.mesh.ny
        bfield = box.sim.interactions[findfirst(x -> isa(x, Demag),
                                                box.sim.interactions)].field
        cnx, cny, cnz = amr.dims[l - 1]
        tx = rx ? (max(1, cld(box.ir[1], 2)):min(cnx, fld(box.ir[end], 2))) : box.ir
        ty = ry ? (max(1, cld(box.jr[1], 2)):min(cny, fld(box.jr[end], 2))) : box.jr
        tz = rz ? (max(1, cld(box.kr[1], 2)):min(cnz, fld(box.kr[end], 2))) : box.kr
        co! = amr_box_to_parent_kernel!(get_backend(amr.Phi[l - 1]), groupsize[])
        co!(amr.Phi[l - 1], bfield, box.ir[1], box.jr[1], box.kr[1], ngx, ngy,
            cnx, cny, first(tx), first(ty), first(tz), rx, ry, rz;
            ndrange=(length(tx), length(ty), length(tz)))
    end
    return nothing
end
