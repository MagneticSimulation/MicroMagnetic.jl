# Flagging and composite-grid generation following Section III of
# García-Cervera & Roma (2006):
#   * flag cells where |div m| >= refine_threshold * max|div m| (the paper's
#     in-plane divergence criterion, eq. (8); with nz = 1 the divergence is
#     in-plane automatically);
#   * optionally flag a layer of cells close to the sample boundary so that
#     the boundary region is covered by the finest level;
#   * cluster the flags into rectangular patches with the Berger-Rigoutsos
#     point-clustering algorithm at a prescribed grid efficiency;
#   * regenerate the composite grid ("remesh") every remesh_interval steps.

"""A rectangle of cells in the coordinates of one level."""
struct AMRRect
    ir::UnitRange{Int}
    jr::UnitRange{Int}
    kr::UnitRange{Int}
end

AMRRect(ir, jr) = AMRRect(ir, jr, 1:1)

_rect_volume(r::AMRRect) = length(r.ir) * length(r.jr) * length(r.kr)

"""Divergence |div m| per cell at every flagging level (1..levels-1),
computed from the composite boxes with central differences (one-sided at the
domain boundary). Returns the arrays and the global maximum."""
function _divergence_arrays(amr::AMRSim{T}) where {T<:AbstractFloat}
    rx, ry, rz = amr.refined
    divs = Vector{Vector{T}}()
    divmax = T(0)
    for l in 1:(amr.levels - 1)
        nx, ny, nz = amr.dims[l]
        hx, hy, hz = _level_h(amr.base_mesh, amr.refined, l)
        C = amr.C[l]
        div = zeros(T, nx * ny * nz)
        for k in 1:nz, j in 1:ny, i in 1:nx
            I = _cell_index(i, j, k, nx, ny)
            d = T(0)
            if rx && nx > 1
                ip = i == nx ? i : i + 1
                im = i == 1 ? i : i - 1
                d += (C[3 * _cell_index(ip, j, k, nx, ny) - 2] -
                      C[3 * _cell_index(im, j, k, nx, ny) - 2]) / ((ip - im) * hx)
            end
            if ry && ny > 1
                jp = j == ny ? j : j + 1
                jm = j == 1 ? j : j - 1
                d += (C[3 * _cell_index(i, jp, k, nx, ny) - 1] -
                      C[3 * _cell_index(i, jm, k, nx, ny) - 1]) / ((jp - jm) * hy)
            end
            if rz && nz > 1
                kp = k == nz ? k : k + 1
                km = k == 1 ? k : k - 1
                d += (C[3 * _cell_index(i, j, kp, nx, ny)] -
                      C[3 * _cell_index(i, j, km, nx, ny)]) / ((kp - km) * hz)
            end
            div[I] = abs(d)
            div[I] > divmax && (divmax = div[I])
        end
        push!(divs, div)
    end
    return divs, divmax
end

"""Flag cells of level `l` (divergence criterion restricted to `valid`, plus
the boundary layer) for the generation of level `l+1` patches."""
function _level_flags(amr::AMRSim{T}, div::Vector{T}, l::Int, theta::T,
                      divmax::T, valid::Vector{Bool}) where {T<:AbstractFloat}
    nx, ny, nz = amr.dims[l]
    flags = zeros(Bool, nx * ny * nz)
    thr = theta * divmax
    if divmax > 0
        for I in eachindex(flags)
            if valid[I] && div[I] >= thr
                flags[I] = true
            end
        end
    end
    if amr.refine_boundary
        bl = amr.boundary_layers
        for k in 1:nz, j in 1:ny, i in 1:nx
            I = _cell_index(i, j, k, nx, ny)
            # outermost layers along every refined dimension
            near_edge = (nx > 1 && (i <= bl || i > nx - bl)) ||
                        (ny > 1 && (j <= bl || j > ny - bl)) ||
                        (nz > 1 && (k <= bl || k > nz - bl))
            if valid[I] && near_edge
                flags[I] = true
            end
        end
    end
    return flags
end

# ------------------------------------------------------- Berger-Rigoutsos

"""Signature of the flags in `box` along dimension `d` (0 = x, 1 = y, 2 = z)."""
function _signature(flags::Vector{Bool}, box::AMRRect, nx::Int, ny::Int, d::Int)
    n = d == 0 ? length(box.ir) : d == 1 ? length(box.jr) : length(box.kr)
    sig = zeros(Int, n)
    for kk in box.kr, jj in box.jr, ii in box.ir
        flags[_cell_index(ii, jj, kk, nx, ny)] || continue
        p = d == 0 ? ii - first(box.ir) + 1 : d == 1 ? jj - first(box.jr) + 1 :
            kk - first(box.kr) + 1
        sig[p] += 1
    end
    return sig
end

"""
Pick a split position in one dimension from its signature, following
Berger & Rigoutsos (1991): prefer the middle of the widest gap (run of zeros),
otherwise a position with maximum slope |dsig| >= 0.95 * max(sig). Returns
`(cut, score)` with `cut = 0` when the signature gives no usable split
(uniform coverage) -- the box is then accepted as is instead of being split
into slivers.
"""
function _sig_cut(sig::Vector{Int})
    n = length(sig)
    n <= 1 && return 0, 0
    M = maximum(sig)
    M == 0 && return 0, 0
    # widest run of zeros
    best_gap, best_gap_cut, run_start = 0, 0, 0
    for p in 1:(n + 1)
        if p <= n && sig[p] == 0
            run_start == 0 && (run_start = p)
        elseif run_start != 0
            len = p - run_start
            if len > best_gap && run_start > 1 && p - 1 < n
                best_gap, best_gap_cut = len, run_start - 1 + (len + 1) >> 1
            end
            run_start = 0
        end
    end
    # maximum slope
    ds, dp = 0, 0
    for p in 1:(n - 1)
        d = abs(sig[p + 1] - sig[p])
        if d > ds
            ds, dp = d, p
        end
    end
    if best_gap_cut > 0
        return best_gap_cut, ds + best_gap
    elseif ds >= 0.95 * M
        return dp, ds
    end
    return 0, 0
end

"""
Cluster the flagged cells inside `box` into rectangles whose flagged fraction
is at least `eff` (grid efficiency). Returns a list of disjoint rectangles
covering all flagged cells.
"""
function _br_cluster(flags::Vector{Bool}, nx::Int, ny::Int, nz::Int,
                     box::AMRRect, eff::Float64)
    total = _rect_volume(box)
    flagged = 0
    for kk in box.kr, jj in box.jr, ii in box.ir
        flags[_cell_index(ii, jj, kk, nx, ny)] && (flagged += 1)
    end
    flagged == 0 && return AMRRect[]
    if flagged / total >= eff || total == 1
        return [box]
    end

    best = (0, 0, 0, -1.0)   # (dim+1, cut, len, score); dim+1 = 0 means no candidate
    for (d, len) in ((0, length(box.ir)), (1, length(box.jr)), (2, length(box.kr)))
        len <= 1 && continue
        sig = _signature(flags, box, nx, ny, d)
        cut, score = _sig_cut(sig)
        cut > 0 && score > best[4] && (best = (d + 1, cut, len, score))
    end
    best[1] == 0 && return [box]

    d, cut = best[1] - 1, best[2]
    if d == 0
        b1 = AMRRect(first(box.ir):(first(box.ir) + cut - 1), box.jr, box.kr)
        b2 = AMRRect((first(box.ir) + cut):last(box.ir), box.jr, box.kr)
    elseif d == 1
        b1 = AMRRect(box.ir, first(box.jr):(first(box.jr) + cut - 1), box.kr)
        b2 = AMRRect(box.ir, (first(box.jr) + cut):last(box.jr), box.kr)
    else
        b1 = AMRRect(box.ir, box.jr, first(box.kr):(first(box.kr) + cut - 1))
        b2 = AMRRect(box.ir, box.jr, (first(box.kr) + cut):last(box.kr))
    end
    out = vcat(_br_cluster(flags, nx, ny, nz, b1, eff),
               _br_cluster(flags, nx, ny, nz, b2, eff))
    return out
end

# ------------------------------------------------------- patch generation

"""Mask of the cells covered by any of `rects` (in the rects' own level)."""
function _rects_mask(rects::Vector{AMRRect}, dims::Tuple{Int,Int,Int})
    mask = zeros(Bool, prod(dims))
    for r in rects, kk in r.kr, jj in r.jr, ii in r.ir
        mask[_cell_index(ii, jj, kk, dims[1], dims[2])] = true
    end
    return mask
end

"""Flag-level rectangles -> level-(l+1) patch rectangles (factor 2 per
refined dimension). Rectangles holding fewer than `min_cells` cells are
dropped: Berger-Rigoutsos recursion can otherwise end in single-cell slivers
whose ghost overhead dominates their useful work."""
function _rects_to_finer(rects::Vector{AMRRect}, refined::NTuple{3,Bool};
                         min_cells::Int=16)
    out = AMRRect[]
    for r in rects
        ir = refined[1] ? ((2 * first(r.ir) - 1):(2 * last(r.ir))) : r.ir
        jr = refined[2] ? ((2 * first(r.jr) - 1):(2 * last(r.jr))) : r.jr
        kr = refined[3] ? ((2 * first(r.kr) - 1):(2 * last(r.kr))) : r.kr
        _rect_volume(AMRRect(ir, jr, kr)) >= min_cells && push!(out, AMRRect(ir, jr, kr))
    end
    return out
end

"""Clip each rectangle in `boxes` into the disjoint parent footprint
`parents` (same-level coordinates): every returned piece lies fully inside
one parent rect. Berger-Rigoutsos boxes are bounding boxes of split regions
and may stick out of the parent patches even when all flagged cells are
inside them; without this clip the next level's patches break proper nesting
and the composite grid double-claims physical cells."""
function _clip_to_parents(boxes::Vector{AMRRect}, parents::Vector{AMRRect})
    out = AMRRect[]
    for b in boxes, p in parents
        ic = intersect(b.ir, p.ir)
        jc = intersect(b.jr, p.jr)
        kc = intersect(b.kr, p.kr)
        (isempty(ic) || isempty(jc) || isempty(kc)) && continue
        push!(out, AMRRect(ic, jc, kc))
    end
    return out
end

"""Generate the rectangles of every patch level (coarse to fine, properly
nested: flags of level l+1 are only taken inside the new level-l patches, and
every clustered box is clipped back into the level-l footprint)."""
function _generate_patches(amr::AMRSim{T}, divs::Vector{Vector{T}}, theta::T,
                           divmax::T) where {T<:AbstractFloat}
    L = amr.levels
    all_rects = Vector{Vector{AMRRect}}()
    valid = ones(Bool, prod(amr.dims[1]))    # level 1 covers the whole domain
    parents = [AMRRect(1:amr.dims[1][1], 1:amr.dims[1][2], 1:amr.dims[1][3])]
    for l in 1:(L - 1)
        flags = _level_flags(amr, divs[l], l, theta, divmax, valid)
        nx, ny, nz = amr.dims[l]
        boxes = any(flags) ?
                _br_cluster(flags, nx, ny, nz,
                            AMRRect(1:nx, 1:ny, 1:nz), Float64(amr.grid_efficiency)) :
                AMRRect[]
        finer = _rects_to_finer(_clip_to_parents(boxes, parents), amr.refined)
        push!(all_rects, finer)
        valid = _rects_mask(finer, amr.dims[l + 1])
        parents = finer
    end
    return all_rects
end

"""
    remesh!(amr)

Regenerate the composite grid from the current magnetization state (which
must be synchronized, i.e. `sync_composite!` run): flag, cluster, rebuild the
patch regions (initialized from the current composite state), rebuild the
GPSM contexts, and refresh masks/ghosts. The divergence threshold is doubled
until the finest-level coverage drops below `max_coverage`.
"""
function remesh!(amr::AMRSim)
    T = eltype(amr.base.sim.spin)
    divs, divmax = _divergence_arrays(amr)
    theta = amr.refine_threshold
    all_rects = Vector{Vector{AMRRect}}()
    for _ in 1:6
        all_rects = _generate_patches(amr, divs, theta, divmax)
        L = amr.levels
        cov = L == 1 ? 0.0 :
              sum(_rect_volume(r) for r in all_rects[L - 1]; init=0) / prod(amr.dims[L])
        cov <= amr.max_coverage && break
        theta *= 2
    end
    if amr.levels > 1
        cov = sum(_rect_volume(r) for r in all_rects[amr.levels - 1]; init=0) /
              prod(amr.dims[amr.levels])
        cov > amr.max_coverage && @warn "AMR remesh: finest-level coverage " *
            "$(round(cov; digits=3)) still exceeds max_coverage=$(amr.max_coverage) " *
            "after doubling the divergence threshold 6 times (theta=$theta); " *
            "accepting the over-covered layout"
    end

    old_auth = [copy(amr.auth[l]) for l in 1:amr.levels]

    # masks from the new rectangles BEFORE the regions are built: the GPSM
    # contexts need the covered masks (Dirichlet cells under finer patches)
    amr.auth[1] .= true
    for l in 2:amr.levels
        amr.auth[l] = _rects_mask(all_rects[l - 1], amr.dims[l])
    end
    for l in 1:(amr.levels - 1)
        fine_mask = amr.auth[l + 1]
        ncx, ncy, ncz = amr.dims[l]
        rx, ry, rz = amr.refined
        nfx, nfy = ncx * (rx ? 2 : 1), ncy * (ry ? 2 : 1)
        cov = amr.covered[l]
        for k in 1:ncz, j in 1:ncy, i in 1:ncx
            ok = true
            for kk in (rz ? 2 * k - 1 : k, rz ? 2 * k : k),
                jj in (ry ? 2 * j - 1 : j, ry ? 2 * j : j),
                ii in (rx ? 2 * i - 1 : i, rx ? 2 * i : i)
                ok &= fine_mask[_cell_index(ii, jj, kk, nfx, nfy)]
            end
            cov[_cell_index(i, j, k, ncx, ncy)] = ok
        end
    end
    amr.covered[amr.levels] .= false

    new_patches = Vector{Vector{AMRRegion{T}}}()
    for l in 2:amr.levels
        ps = AMRRegion{T}[]
        for r in all_rects[l - 1]
            push!(ps, _make_patch(amr, l, r.ir, r.jr, r.kr, old_auth[l]))
        end
        push!(new_patches, ps)
    end
    amr.patches = new_patches

    sync_composite!(amr)
    fill_ghosts!(amr)
    rebuild_boxes!(amr)
    amr.n_remesh += 1
    amr.g_initialized = false   # new regions need a fresh GPSM initialization
    amr.ecache = nothing        # state and geometry changed: energy cache stale
    return nothing
end
