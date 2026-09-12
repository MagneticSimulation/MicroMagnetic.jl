# Composite-grid synchronization: the per-level composite magnetization boxes
# C[l], the composite auxiliary-field boxes Gc[l] of the GPSM integrator, the
# fine->coarse restriction, the coarse->fine interpolation and the ghost-cell
# filling of every patch.
#
# Layout conventions (matching FDMesh): cell (i,j,k) sits at flat index
# (k-1)*nx*ny + (j-1)*nx + i, and the magnetization at 3*I-2 .. 3*I. Patch
# boxes carry ghost layers on the refined dimensions only: box cell (a,b,c)
# with gx < a <= ngx-gx (per dimension, gx = ghost on refined dims, 0
# otherwise) is interior, and interior cell (a-gx, b-gy, c-gz) maps to the
# global level-`level` cell (i0-1+a, j0-1+b, k0-1+c).

# Parent cells and interpolation weights of the cell-centered coarse->fine
# linear interpolation. Coarse cell i has its center at (i-0.5)h, a fine cell
# j at (j-0.5)h/2: odd fine cells lean on the coarse cells (j>>1, j>>1 + 1)
# with weights (1/4, 3/4), even fine cells on (j>>1, j>>1 + 1) with weights
# (3/4, 1/4). At the domain edges the out-of-range parent is folded onto the
# boundary cell (constant extension); unrefined dimensions copy through.
@inline function _pw(rx::Bool, i::Int, np::Int)
    if !rx
        return i, i, 1.0, 0.0
    end
    ic_raw = i >> 1
    ic = clamp(ic_raw, 1, np)
    ic2 = clamp(ic_raw + 1, 1, np)
    w1 = isodd(i) ? 0.25 : 0.75
    return ic, ic2, w1, 1.0 - w1
end

@inline function _gather3(src, i1, i2, w1x, w2x, j1, j2, w1y, w2y, k1, k2, w1z, w2z,
                          npx, npy, c)
    return (w1x * w1y * w1z * src[3 * _cell_index(i1, j1, k1, npx, npy) - 2 + c] +
            w2x * w1y * w1z * src[3 * _cell_index(i2, j1, k1, npx, npy) - 2 + c] +
            w1x * w2y * w1z * src[3 * _cell_index(i1, j2, k1, npx, npy) - 2 + c] +
            w2x * w2y * w1z * src[3 * _cell_index(i2, j2, k1, npx, npy) - 2 + c] +
            w1x * w1y * w2z * src[3 * _cell_index(i1, j1, k2, npx, npy) - 2 + c] +
            w2x * w1y * w2z * src[3 * _cell_index(i2, j1, k2, npx, npy) - 2 + c] +
            w1x * w2y * w2z * src[3 * _cell_index(i1, j2, k2, npx, npy) - 2 + c] +
            w2x * w2y * w2z * src[3 * _cell_index(i2, j2, k2, npx, npy) - 2 + c])
end

"""
Fill the whole level-`ell` box `dst` with the coarse->fine interpolation of
the level-`ell-1` box `src`. Both boxes are globally aligned. With
`normalize = true` the interpolated magnetization (a direction field) is
renormalized; the GPSM auxiliary field g is transferred with
`normalize = false`.
"""
@kernel function amr_interp_box_kernel!(dst, @Const(src), rx::Bool, ry::Bool, rz::Bool,
                                        ndx::Int, ndy::Int, ndz::Int,
                                        npx::Int, npy::Int, npz::Int, normalize::Bool)
    i, j, k = @index(Global, NTuple)
    i1, i2, wx1, wx2 = _pw(rx, i, npx)
    j1, j2, wy1, wy2 = _pw(ry, j, npy)
    k1, k2, wz1, wz2 = _pw(rz, k, npz)
    d = 3 * _cell_index(i, j, k, ndx, ndy) - 2
    @inbounds begin
        mx = _gather3(src, i1, i2, wx1, wx2, j1, j2, wy1, wy2, k1, k2, wz1, wz2,
                      npx, npy, 0)
        my = _gather3(src, i1, i2, wx1, wx2, j1, j2, wy1, wy2, k1, k2, wz1, wz2,
                      npx, npy, 1)
        mz = _gather3(src, i1, i2, wx1, wx2, j1, j2, wy1, wy2, k1, k2, wz1, wz2,
                      npx, npy, 2)
        if normalize
            nm = sqrt(mx * mx + my * my + mz * mz)
            if nm > 1e-8
                mx /= nm
                my /= nm
                mz /= nm
            end
        end
        dst[d] = mx
        dst[d + 1] = my
        dst[d + 2] = mz
    end
end

"""
Fill the ghost ring of a patch box by interpolating the coarser composite box
`src` (level `level-1`) at the ghost-cell global positions. Cells outside the
sample are filled by constant extension of the boundary values, which makes
the exchange stencil across the sample boundary a no-op (free boundary,
matching the ngbs = -1 treatment of uniform sims). The interpolated
magnetization is renormalized.
"""
@kernel function amr_fill_ghosts_kernel!(box, @Const(src), rx::Bool, ry::Bool, rz::Bool,
                                         gx::Int, gy::Int, gz::Int,
                                         i0::Int, j0::Int, k0::Int,
                                         ngx::Int, ngy::Int, ngz::Int,
                                         nlx::Int, nly::Int, nlz::Int,
                                         npx::Int, npy::Int, npz::Int)
    a, b, c = @index(Global, NTuple)
    interior = gx < a <= ngx - gx && gy < b <= ngy - gy && gz < c <= ngz - gz
    if !interior
        gi = clamp(i0 - 1 + a, 1, nlx)
        gj = clamp(j0 - 1 + b, 1, nly)
        gk = clamp(k0 - 1 + c, 1, nlz)
        i1, i2, wx1, wx2 = _pw(rx, gi, npx)
        j1, j2, wy1, wy2 = _pw(ry, gj, npy)
        k1, k2, wz1, wz2 = _pw(rz, gk, npz)
        d = 3 * _cell_index(a, b, c, ngx, ngy) - 2
        @inbounds begin
            mx = _gather3(src, i1, i2, wx1, wx2, j1, j2, wy1, wy2, k1, k2, wz1, wz2,
                          npx, npy, 0)
            my = _gather3(src, i1, i2, wx1, wx2, j1, j2, wy1, wy2, k1, k2, wz1, wz2,
                          npx, npy, 1)
            mz = _gather3(src, i1, i2, wx1, wx2, j1, j2, wy1, wy2, k1, k2, wz1, wz2,
                          npx, npy, 2)
            nm = sqrt(mx * mx + my * my + mz * mz)
            if nm > 1e-8
                mx /= nm
                my /= nm
                mz /= nm
            end
            box[d] = mx
            box[d + 1] = my
            box[d + 2] = mz
        end
    end
end

"""Copy the interior of a patch box into the composite box (patch -> C)."""
@kernel function amr_patch_to_box_kernel!(box, @Const(patch), gx::Int, gy::Int, gz::Int,
                                          i0::Int, j0::Int,
                                          k0::Int, ngx::Int, ngy::Int,
                                          nbx::Int, nby::Int)
    a, b, c = @index(Global, NTuple)
    Ib = _cell_index(a + gx, b + gy, c + gz, ngx, ngy)
    Ig = _cell_index(i0 - 1 + a, j0 - 1 + b, k0 - 1 + c, nbx, nby)
    d = 3 * Ib - 2
    e = 3 * Ig - 2
    @inbounds begin
        box[e] = patch[d]
        box[e + 1] = patch[d + 1]
        box[e + 2] = patch[d + 2]
    end
end

"""Copy the composite box into the interior of a patch box (C -> patch)."""
@kernel function amr_box_to_patch_kernel!(patch, @Const(box), gx::Int, gy::Int, gz::Int,
                                          i0::Int, j0::Int,
                                          k0::Int, ngx::Int, ngy::Int,
                                          nbx::Int, nby::Int)
    a, b, c = @index(Global, NTuple)
    Ib = _cell_index(a + gx, b + gy, c + gz, ngx, ngy)
    Ig = _cell_index(i0 - 1 + a, j0 - 1 + b, k0 - 1 + c, nbx, nby)
    d = 3 * Ib - 2
    e = 3 * Ig - 2
    @inbounds begin
        patch[d] = box[e]
        patch[d + 1] = box[e + 1]
        patch[d + 2] = box[e + 2]
    end
end

"""
Restrict the level `ell+1` box onto the cells of level `ell` that lie under
finer patches (`covered`): each covered coarse cell receives the average of
its (2^r) children. With `normalize = true` (magnetization) the average is
renormalized to unit length; the GPSM auxiliary field uses `normalize = false`.
"""
@kernel function amr_restrict_kernel!(coarse, @Const(fine), rx::Bool, ry::Bool, rz::Bool,
                                      ncx::Int, ncy::Int, @Const(covered), normalize::Bool)
    i, j, k = @index(Global, NTuple)
    I = _cell_index(i, j, k, ncx, ncy)
    if covered[I]
        i1 = rx ? 2 * i - 1 : i
        i2 = rx ? 2 * i : i
        j1 = ry ? 2 * j - 1 : j
        j2 = ry ? 2 * j : j
        k1 = rz ? 2 * k - 1 : k
        k2 = rz ? 2 * k : k
        n = 0.0
        mx = 0.0
        my = 0.0
        mz = 0.0
        @inbounds for kk in (k1, k2), jj in (j1, j2), ii in (i1, i2)
            e = 3 * _cell_index(ii, jj, kk, ncx * (rx ? 2 : 1), ncy * (ry ? 2 : 1)) - 2
            mx += fine[e]
            my += fine[e + 1]
            mz += fine[e + 2]
            n += 1.0
        end
        mx /= n
        my /= n
        mz /= n
        if normalize
            nm = sqrt(mx * mx + my * my + mz * mz)
            if nm > 1e-8
                mx /= nm
                my /= nm
                mz /= nm
            end
        end
        d = 3 * I - 2
        @inbounds begin
            coarse[d] = mx
            coarse[d + 1] = my
            coarse[d + 2] = mz
        end
    end
end

"""
Initialize a new patch interior from the current composite state: cells that
were already authoritative at this level (old fine data, per `fine_auth`) are
copied from the level-`l` box, the rest are interpolated from level `l-1`
(renormalized).
"""
@kernel function amr_init_patch_kernel!(box, @Const(fine_box), @Const(fine_auth),
                                        @Const(coarse_box), use_fine::Bool,
                                        rx::Bool, ry::Bool, rz::Bool,
                                        gx::Int, gy::Int, gz::Int,
                                        i0::Int, j0::Int, k0::Int,
                                        ngx::Int, ngy::Int,
                                        nlx::Int, nly::Int, nlz::Int,
                                        npx::Int, npy::Int, npz::Int)
    a, b, c = @index(Global, NTuple)
    gi = clamp(i0 - 1 + a, 1, nlx)
    gj = clamp(j0 - 1 + b, 1, nly)
    gk = clamp(k0 - 1 + c, 1, nlz)
    Ib = _cell_index(a + gx, b + gy, c + gz, ngx, ngy)
    Ig = _cell_index(gi, gj, gk, nlx, nly)
    d = 3 * Ib - 2
    e = 3 * Ig - 2
    if use_fine && fine_auth[Ig]
        @inbounds begin
            box[d] = fine_box[e]
            box[d + 1] = fine_box[e + 1]
            box[d + 2] = fine_box[e + 2]
        end
    else
        i1, i2, wx1, wx2 = _pw(rx, gi, npx)
        j1, j2, wy1, wy2 = _pw(ry, gj, npy)
        k1, k2, wz1, wz2 = _pw(rz, gk, npz)
        @inbounds begin
            mx = _gather3(coarse_box, i1, i2, wx1, wx2, j1, j2, wy1, wy2,
                          k1, k2, wz1, wz2, npx, npy, 0)
            my = _gather3(coarse_box, i1, i2, wx1, wx2, j1, j2, wy1, wy2,
                          k1, k2, wz1, wz2, npx, npy, 1)
            mz = _gather3(coarse_box, i1, i2, wx1, wx2, j1, j2, wy1, wy2,
                          k1, k2, wz1, wz2, npx, npy, 2)
            nm = sqrt(mx * mx + my * my + mz * mz)
            if nm > 1e-8
                mx /= nm
                my /= nm
                mz /= nm
            end
            box[d] = mx
            box[d + 1] = my
            box[d + 2] = mz
        end
    end
end

"""Scatter the three GPSM auxiliary component arrays of a region into the
3N-interleaved composite box (g1,g2,g3 -> Gc)."""
@kernel function amr_g_to_box_kernel!(box, @Const(g1), @Const(g2), @Const(g3),
                                      gx::Int, gy::Int, gz::Int, i0::Int, j0::Int,
                                      k0::Int, ngx::Int, ngy::Int, nbx::Int, nby::Int)
    a, b, c = @index(Global, NTuple)
    Ib = _cell_index(a + gx, b + gy, c + gz, ngx, ngy)
    Ig = _cell_index(i0 - 1 + a, j0 - 1 + b, k0 - 1 + c, nbx, nby)
    d = 3 * Ig - 2
    @inbounds begin
        box[d] = g1[Ib]
        box[d + 1] = g2[Ib]
        box[d + 2] = g3[Ib]
    end
end

"""Gather the 3N-interleaved composite auxiliary box into a region's three
component arrays (Gc -> g1,g2,g3)."""
@kernel function amr_box_to_g_kernel!(g1, g2, g3, @Const(box), gx::Int, gy::Int, gz::Int,
                                      i0::Int, j0::Int,
                                      k0::Int, ngx::Int, ngy::Int, nbx::Int, nby::Int)
    a, b, c = @index(Global, NTuple)
    Ib = _cell_index(a + gx, b + gy, c + gz, ngx, ngy)
    Ig = _cell_index(i0 - 1 + a, j0 - 1 + b, k0 - 1 + c, nbx, nby)
    d = 3 * Ib - 2
    e = 3 * Ig - 2
    @inbounds begin
        g1[Ib] = box[e]
        g2[Ib] = box[e + 1]
        g3[Ib] = box[e + 2]
    end
end

# ---------------------------------------------------------------- drivers

"""Ghost-layer count per dimension: only refined dimensions carry ghosts."""
_ghosts(amr::AMRSim) = (amr.refined[1] ? amr.ghost : 0,
                        amr.refined[2] ? amr.ghost : 0,
                        amr.refined[3] ? amr.ghost : 0)

"""The interior cell count of a patch box (host side)."""
_patch_interior_dims(p::AMRRegion) =
    (length(p.ir), length(p.jr), length(p.kr))

function _up_pass!(amr::AMRSim)
    rx, ry, rz = amr.refined
    for l in 2:amr.levels
        ndx, ndy, ndz = amr.dims[l]
        npx, npy, npz = amr.dims[l - 1]
        kernel! = amr_interp_box_kernel!(get_backend(amr.C[l]), groupsize[])
        kernel!(amr.C[l], amr.C[l - 1], rx, ry, rz, ndx, ndy, ndz, npx, npy, npz, true;
                ndrange=(ndx, ndy, ndz))
        for p in amr.patches[l - 1]
            g = p.geo
            ck! = amr_patch_to_box_kernel!(get_backend(amr.C[l]), groupsize[])
            ck!(amr.C[l], p.sim.spin, g.gx, g.gy, g.gz,
                g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
                g.ngx, g.ngy, g.nlx, g.nly;
                ndrange=(g.nxi, g.nyi, g.nzi))
        end
    end
end

function _down_pass!(amr::AMRSim)
    rx, ry, rz = amr.refined
    for l in (amr.levels - 1):-1:1
        ncx, ncy, ncz = amr.dims[l]
        kernel! = amr_restrict_kernel!(get_backend(amr.C[l]), groupsize[])
        kernel!(amr.C[l], amr.C[l + 1], rx, ry, rz, ncx, ncy, amr.covered[l], true;
                ndrange=(ncx, ncy, ncz))
        if l >= 2
            for p in amr.patches[l - 1]
                g = p.geo
                bk! = amr_box_to_patch_kernel!(get_backend(p.sim.spin), groupsize[])
                bk!(p.sim.spin, amr.C[l], g.gx, g.gy, g.gz,
                    g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
                    g.ngx, g.ngy, g.nlx, g.nly;
                    ndrange=(g.nxi, g.nyi, g.nzi))
            end
        end
    end
end

"""
    sync_composite!(amr)

Bring the per-level composite magnetization boxes `amr.C` in sync with the
authoritative region data (base spin + patch interiors): an up pass builds
every level from the coarser one and overwrites it with patch data, a down
pass averages fine data onto the coarse cells under finer patches (writing
the restriction back into the owning patches and into the base spin), and a
second up pass refreshes the interpolated regions.
"""
function sync_composite!(amr::AMRSim)
    _up_pass!(amr)
    _down_pass!(amr)
    _up_pass!(amr)
    return nothing
end

"""
    sync_g_boxes!(amr)

Same three-pass synchronization for the composite GPSM auxiliary field boxes
`amr.Gc`, WITHOUT renormalization (g is not a direction field). Before the
first solve (`!g_initialized`) the auxiliary field is approximated by the
composite magnetization.
"""
function sync_g_boxes!(amr::AMRSim)
    if !amr.g_initialized
        for l in 1:amr.levels
            copyto!(amr.Gc[l], amr.C[l])
        end
        return nothing
    end
    rx, ry, rz = amr.refined
    # up: Gc[l] = P(Gc[l-1]) + patch g
    for l in 2:amr.levels
        ndx, ndy, ndz = amr.dims[l]
        npx, npy, npz = amr.dims[l - 1]
        kernel! = amr_interp_box_kernel!(get_backend(amr.Gc[l]), groupsize[])
        kernel!(amr.Gc[l], amr.Gc[l - 1], rx, ry, rz, ndx, ndy, ndz, npx, npy, npz,
                false; ndrange=(ndx, ndy, ndz))
        for p in amr.patches[l - 1]
            g = p.geo
            gk! = amr_g_to_box_kernel!(get_backend(amr.Gc[l]), groupsize[])
            gk!(amr.Gc[l], p.g1, p.g2, p.g3, g.gx, g.gy, g.gz,
                g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
                g.ngx, g.ngy, g.nlx, g.nly;
                ndrange=(g.nxi, g.nyi, g.nzi))
        end
    end
    # down: restricted fine g onto covered coarse cells (+ back into patches)
    for l in (amr.levels - 1):-1:1
        ncx, ncy, ncz = amr.dims[l]
        if l == 1
            # the base region's own g enters the composite box first
            bg = amr.base.geo
            bx, by, bz = amr.dims[1]
            gk! = amr_g_to_box_kernel!(get_backend(amr.Gc[1]), groupsize[])
            gk!(amr.Gc[1], amr.base.g1, amr.base.g2, amr.base.g3,
                bg.gx, bg.gy, bg.gz, bg.i0 + bg.gx, bg.j0 + bg.gy, bg.k0 + bg.gz,
                bg.ngx, bg.ngy, bg.nlx, bg.nly; ndrange=(bx, by, bz))
        end
        kernel! = amr_restrict_kernel!(get_backend(amr.Gc[l]), groupsize[])
        kernel!(amr.Gc[l], amr.Gc[l + 1], rx, ry, rz, ncx, ncy, amr.covered[l], false;
                ndrange=(ncx, ncy, ncz))
        if l >= 2
            for p in amr.patches[l - 1]
                g = p.geo
                bg! = amr_box_to_g_kernel!(get_backend(p.g1), groupsize[])
                bg!(p.g1, p.g2, p.g3, amr.Gc[l], g.gx, g.gy, g.gz,
                    g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
                    g.ngx, g.ngy, g.nlx, g.nly;
                    ndrange=(g.nxi, g.nyi, g.nzi))
            end
        else
            # level 1: the base region's own g arrays follow the composite box
            bg = amr.base.geo
            bx, by, bz = amr.dims[1]
            bg! = amr_box_to_g_kernel!(get_backend(amr.base.g1), groupsize[])
            bg!(amr.base.g1, amr.base.g2, amr.base.g3, amr.Gc[1],
                bg.gx, bg.gy, bg.gz, bg.i0 + bg.gx, bg.j0 + bg.gy, bg.k0 + bg.gz,
                bg.ngx, bg.ngy, bg.nlx, bg.nly; ndrange=(bx, by, bz))
        end
    end
    # up again with the updated coarser data
    for l in 2:amr.levels
        ndx, ndy, ndz = amr.dims[l]
        npx, npy, npz = amr.dims[l - 1]
        kernel! = amr_interp_box_kernel!(get_backend(amr.Gc[l]), groupsize[])
        kernel!(amr.Gc[l], amr.Gc[l - 1], rx, ry, rz, ndx, ndy, ndz, npx, npy, npz,
                false; ndrange=(ndx, ndy, ndz))
        for p in amr.patches[l - 1]
            g = p.geo
            gk! = amr_g_to_box_kernel!(get_backend(amr.Gc[l]), groupsize[])
            gk!(amr.Gc[l], p.g1, p.g2, p.g3, g.gx, g.gy, g.gz,
                g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
                g.ngx, g.ngy, g.nlx, g.nly;
                ndrange=(g.nxi, g.nyi, g.nzi))
        end
    end
    return nothing
end

"""
    fill_ghosts!(amr)

Refresh the ghost cells of every patch from the next coarser composite box.
"""
function fill_ghosts!(amr::AMRSim)
    rx, ry, rz = amr.refined
    for l in 2:amr.levels
        npx, npy, npz = amr.dims[l - 1]
        for p in amr.patches[l - 1]
            g = p.geo
            kernel! = amr_fill_ghosts_kernel!(get_backend(p.sim.spin), groupsize[])
            kernel!(p.sim.spin, amr.C[l - 1], rx, ry, rz,
                    g.gx, g.gy, g.gz,
                    g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
                    g.ngx, g.ngy, g.ngz, g.nlx, g.nly, g.nlz, npx, npy, npz;
                    ndrange=(g.ngx, g.ngy, g.ngz))
        end
    end
    return nothing
end

"""Initialize a patch interior from the current composite boxes: cells that
were already authoritative at this level (old fine data, per `old_auth`) are
copied from the level-`l` box, the rest are interpolated from level `l-1`."""
function _init_patch_interior!(amr::AMRSim{T}, r::AMRRegion{T},
                               old_auth::Vector{Bool}) where {T<:AbstractFloat}
    rx, ry, rz = amr.refined
    l = r.level
    npx, npy, npz = amr.dims[l - 1]
    use_fine = l <= length(amr.auth) && any(old_auth)
    g = r.geo
    kernel! = amr_init_patch_kernel!(get_backend(r.sim.spin), groupsize[])
    kernel!(r.sim.spin, amr.C[l], old_auth, amr.C[l - 1], use_fine, rx, ry, rz,
            g.gx, g.gy, g.gz,
            g.i0 + g.gx, g.j0 + g.gy, g.k0 + g.gz,
            g.ngx, g.ngy, g.nlx, g.nly, g.nlz, npx, npy, npz;
            ndrange=(g.nxi, g.nyi, g.nzi))
    return nothing
end
