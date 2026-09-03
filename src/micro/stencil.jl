# Stencil acceleration layer for Exchange/DMI.
#
# When the exchange/DMI parameters are per-region uniform (the common case for
# materials set up via set_region + region-mapped constants, or single-material
# samples), the harmonic-mean coupling 2*a*b/(a+b) takes only (R+1)^2 distinct
# values (R = #material classes; class 0 = vacuum). We precompute those into a
# pair table and let a "partition" kernel look it up by (class_I, class_nb)
# instead of gathering per-cell A/D and dividing each bond.
#
# Design notes:
#  - inline kernels (kernels.jl exchange_kernel!/bulkdmi_kernel!/interfacial_dmi_kernel!)
#    are NOT modified. The partition path is a separate fast path taken only when
#    the per-class uniformity scan succeeds; otherwise we fall back to inline.
#  - The pair table is (R+1)x(R+1); class 0 (vacuum) row/column is identically
#    zero, so a vacuum neighbour auto-contributes 0 with no per-bond branch.
#  - Tables are built in Float64 on the host, then converted to T. The <=1 ulp
#    difference vs inline is acceptable (we do not require bit-identity).
#  - Caches key on mesh.layout_version (bumped by set_region) and are dropped
#    directly by set_Ms (Ms redefines the vacuum mask, hence the class map).

using KernelAbstractions

# Keep pair tables within L1 cache; beyond this we fall back to inline rather
# than thrash cache with a large table.
const _MAX_CLASSES = 16

# --------------------------------------------------------------------------------------
# Cache invalidation
# --------------------------------------------------------------------------------------

"""
    _drop_ms_caches!(sim::MicroSim)

Drop all Ms-dependent stencil caches: the class map itself plus every
interaction's pair table. Called by set_Ms because a new Ms redefines which
cells are vacuum, so the class map and all derived tables are stale.

Layout-version-keyed caches (only region changes) are handled per-interaction
in _exchange_tables!/_dmi_tables!, not here.
"""
function _drop_ms_caches!(sim::MicroSim)
    sim.mat_class = nothing
    sim.mat_class_layout = -1
    sim.n_classes = 0
    for E in sim.interactions
        if isa(E, Exchange)
            E.pair_Ax = E.pair_Ay = E.pair_Az = nothing
            E.stencil_layout = -1
        elseif isa(E, DMI)
            E.pair_Dx = E.pair_Dy = E.pair_Dz = nothing
            empty!(E.Dcls)
            E.stencil_layout = -1
        end
    end
    return nothing
end

"""
    _update_inv_ms!(sim::MicroSim)

Recompute sim.inv_ms (1/mu0_Ms, with 0 for vacuum) from sim.mu0_Ms. Stored in
precision T so partition kernels divide by multiplying the precomputed inverse
(no per-cell division, no Float64 promotion). Called by set_Ms.
"""
function _update_inv_ms!(sim::MicroSim)
    T = eltype(sim.mu0_Ms)
    ms = sim.mu0_Ms
    if ms isa Fill
        v = ms.value
        sim.inv_ms = Fill(v > 0 ? inv(T(v)) : T(0), length(ms))
    else
        # build on host in T directly (inv fits in T for typical Ms magnitudes);
        # 0 stays 0 to mark vacuum.
        a = Array(ms)
        out = Array{T,1}(undef, length(a))
        @inbounds @simd for i in eachindex(a)
            v = a[i]
            out[i] = v > 0 ? inv(T(v)) : T(0)
        end
        sim.inv_ms = kernel_array(out)
    end
    return nothing
end

# --------------------------------------------------------------------------------------
# Class map: region + Ms -> compact material classes (0 = vacuum, 1..R = materials)
# --------------------------------------------------------------------------------------

"""
    _fresh_mat_class(sim) -> (mat_class, n_classes, R_vacuum_ok)

Build (or return cached) per-cell material class map:
  - 0  = vacuum (mu0_Ms == 0), regardless of region
  - 1..R = first-encounter-compressed region id among non-vacuum cells

Class 0 lives in the map so the pair table's class-0 row/column (identically
zero) can mask vacuum bonds with no per-bond branch. Returns:
  - mat_class: Fill or dense UInt8 array (backend kernel array if dense)
  - n_classes: R (number of material classes)
  - R_vacuum_ok: whether the dense (R<=_MAX_CLASSES) path is usable for the
    partition fast path. When R > _MAX_CLASSES the map degrades to a 0/1 vacuum
    mask (mat_class in {0,1}, n_classes=1) and callers must fall back to inline.
"""
function _fresh_mat_class(sim::MicroSim)
    mesh = sim.mesh
    lv = mesh.layout_version[]
    if sim.mat_class !== nothing && sim.mat_class_layout == lv
        cls = sim.mat_class
        if cls isa Fill
            # Fill(0) => all vacuum; Fill(1) => single material, no vacuum.
            return (cls, sim.n_classes, sim.n_classes <= _MAX_CLASSES)
        end
        return (cls, sim.n_classes, sim.n_classes <= _MAX_CLASSES)
    end

    n = sim.n_total
    ms = sim.mu0_Ms
    regions = mesh.regions
    if ms isa Fill
        if ms.value == 0
            # all vacuum: class 0 everywhere
            sim.mat_class = Fill(UInt8(0), n)
            sim.n_classes = 0
            sim.mat_class_layout = lv
            return (sim.mat_class, 0, true)
        end
        # non-vacuum everywhere: classes collapse to regions; but if regions are
        # all 0 (no set_region called), it's a single material -> class 1.
        r = Array(regions)
        distinct = Int[]
        @inbounds for i in 1:n
            ri = r[i]
            pushfirst_if_new!(distinct, ri)
        end
        R = length(distinct)
        if R == 1
            sim.mat_class = Fill(UInt8(1), n)
            sim.n_classes = 1
        elseif R <= _MAX_CLASSES
            # map region -> 1..R by first-encounter order
            map = Dict{Int,UInt8}()
            for (i, rid) in enumerate(distinct)
                map[rid] = UInt8(i)
            end
            cls = Array{UInt8,1}(undef, n)
            @inbounds @simd for i in 1:n
                cls[i] = map[r[i]]
            end
            sim.mat_class = kernel_array(cls)
            sim.n_classes = R
        else
            # too many regions; degrade to vacuum mask (R=1 is a lie, but the
            # partition path will be rejected by the R<=_MAX_CLASSES check).
            sim.mat_class = Fill(UInt8(1), n)
            sim.n_classes = R  # truthfully report R; caller rejects via R_ok
        end
        sim.mat_class_layout = lv
        return (sim.mat_class, sim.n_classes, sim.n_classes <= _MAX_CLASSES)
    end

    # general path: vacuum mask from Ms + region compression among non-vacuum
    a = Array(ms)
    r = Array(regions)
    distinct = Int[]
    @inbounds for i in 1:n
        a[i] == 0 && continue
        pushfirst_if_new!(distinct, r[i])
    end
    R = length(distinct)
    if R == 0
        sim.mat_class = Fill(UInt8(0), n)
        sim.n_classes = 0
        sim.mat_class_layout = lv
        return (sim.mat_class, 0, true)
    elseif R == 1
        # single material; vacuum cells still 0, material cells 1
        cls = Array{UInt8,1}(undef, n)
        @inbounds @simd for i in 1:n
            cls[i] = a[i] > 0 ? UInt8(1) : UInt8(0)
        end
        sim.mat_class = kernel_array(cls)
        sim.n_classes = 1
        sim.mat_class_layout = lv
        return (sim.mat_class, 1, true)
    elseif R <= _MAX_CLASSES
        map = Dict{Int,UInt8}()
        for (i, rid) in enumerate(distinct)
            map[rid] = UInt8(i)
        end
        cls = Array{UInt8,1}(undef, n)
        @inbounds @simd for i in 1:n
            cls[i] = a[i] > 0 ? map[r[i]] : UInt8(0)
        end
        sim.mat_class = kernel_array(cls)
        sim.n_classes = R
        sim.mat_class_layout = lv
        return (sim.mat_class, R, true)
    else
        # degrade to 0/1 vacuum mask; partition path rejected by R_ok
        cls = Array{UInt8,1}(undef, n)
        @inbounds @simd for i in 1:n
            cls[i] = a[i] > 0 ? UInt8(1) : UInt8(0)
        end
        sim.mat_class = kernel_array(cls)
        sim.n_classes = R
        sim.mat_class_layout = lv
        return (sim.mat_class, R, false)
    end
end

# helper: prepend rid to distinct if not already present, keeping first-encounter order
@inline function pushfirst_if_new!(distinct::Vector{Int}, rid::Integer)
    @inbounds for k in eachindex(distinct)
        distinct[k] == rid && return nothing
    end
    pushfirst!(distinct, Int(rid))
    return nothing
end

# --------------------------------------------------------------------------------------
# Per-class uniformity scan
# --------------------------------------------------------------------------------------

"""
    _class_uniform_values(sim, a) -> Union{Nothing,Vector{T}}

If `a` (a parameter array per cell) is constant within each non-vacuum class,
return a length-(R+1) Vector{T} with index (class+1): class 0 = 0, classes 1..R
their uniform value. Otherwise return nothing.

Fill inputs are trivially uniform across all classes and short-circuit.
"""
function _class_uniform_values(sim::MicroSim, a::AbstractArray)
    T = eltype(sim.spin)
    cls, R, R_ok = _fresh_mat_class(sim)
    R_ok || return nothing
    if a isa Fill
        # uniform across everything
        vals = zeros(T, R + 1)
        v = T(a.value)
        @inbounds for c in 1:R
            vals[c+1] = v
        end
        return vals
    end

    # dense: scan once on host; bail on first within-class mismatch
    a_h = Array(a)
    if cls isa Fill
        c0 = cls.value  # 0 => all vacuum (R=0), 1 => single material
        @inbounds if c0 == 0
            return zeros(T, 1)  # all vacuum; trivially uniform
        else
            v0 = a_h[1]
            @inbounds @simd for i in eachindex(a_h)
                if a_h[i] != v0
                    return nothing
                end
            end
            vals = zeros(T, 2)
            vals[2] = T(v0)
            return vals
        end
    end

    cls_h = Array(cls)
    vals = zeros(T, R + 1)
    seen = falses(R + 1)  # seen[c+1]
    @inbounds for i in eachindex(a_h)
        c = cls_h[i]
        c == 0 && continue  # vacuum
        v = a_h[i]
        idx = c + 1
        if !seen[idx]
            seen[idx] = true
            vals[idx] = T(v)
        elseif vals[idx] != T(v)
            return nothing
        end
    end
    return vals
end

# --------------------------------------------------------------------------------------
# Pair table construction (Float64 on host -> T)
# --------------------------------------------------------------------------------------

"""
    _pair_table(vals) -> Matrix{T} (size (R+1) x (R+1))

Harmonic-mean coupling table: pair[cI+1, cNb+1] = 2*a*b/(a+b), 0 if either is 0.
Class 0 (vacuum) row/column is identically zero so vacuum bonds auto-vanish.
Built in Float64 then converted to T to keep the <=1 ulp inline-vs-partition
difference small.
"""
function _pair_table(vals::Vector{T}) where {T<:AbstractFloat}
    R = length(vals) - 1
    pair = zeros(T, R + 1, R + 1)
    @inbounds for cI in 1:R
        a = Float64(vals[cI+1])
        a == 0 && continue
        for cNb in 1:R
            b = Float64(vals[cNb+1])
            b == 0 && continue
            pair[cI+1, cNb+1] = T(2 * a * b / (a + b))
        end
    end
    return pair
end

# --------------------------------------------------------------------------------------
# Interaction table builders (lazy + layout-stamped)
# --------------------------------------------------------------------------------------

"""
    _exchange_tables!(sim, exch) -> (mat_class, ok::Bool)

Ensure exch.pair_A{y,z} are fresh for the current mesh layout. Returns the class
map and whether the partition fast path is usable. On success the three pair
tables are populated (in T, on the chosen backend). On failure they are cleared.
"""
function _exchange_tables!(sim::MicroSim, exch::Exchange{T}) where {T<:AbstractFloat}
    mesh = sim.mesh
    lv = mesh.layout_version[]
    if exch.pair_Ax !== nothing && exch.stencil_layout == lv
        cls, R, R_ok = _fresh_mat_class(sim)
        return (cls, R_ok)
    end

    cls, R, R_ok = _fresh_mat_class(sim)
    R_ok || return (cls, false)

    vx = _class_uniform_values(sim, exch.Ax)
    if vx === nothing
        exch.pair_Ax = exch.pair_Ay = exch.pair_Az = nothing
        exch.stencil_layout = -1
        return (cls, false)
    end
    vy = _class_uniform_values(sim, exch.Ay)
    if vy === nothing
        exch.pair_Ax = exch.pair_Ay = exch.pair_Az = nothing
        exch.stencil_layout = -1
        return (cls, false)
    end
    vz = _class_uniform_values(sim, exch.Az)
    if vz === nothing
        exch.pair_Ax = exch.pair_Ay = exch.pair_Az = nothing
        exch.stencil_layout = -1
        return (cls, false)
    end

    # build pair tables in Float64 -> T, then push to backend
    pair_x = kernel_array(_pair_table(vx))
    pair_y = kernel_array(_pair_table(vy))
    pair_z = kernel_array(_pair_table(vz))
    exch.pair_Ax = pair_x
    exch.pair_Ay = pair_y
    exch.pair_Az = pair_z
    exch.stencil_layout = lv
    return (cls, true)
end

"""
    _dmi_tables!(sim, dmi) -> (mat_class, ok::Bool)

Ensure dmi.pair_D{y,z} (and Dcls for interfacial) are fresh. Time-modulated DMI
(ft !== _static_time) always returns false: the pair table is a static value and
cannot be scaled by tfac.
"""
function _dmi_tables!(sim::MicroSim, dmi::DMI{T}) where {T<:AbstractFloat}
    # time-modulated DMI cannot use a static pair table
    if dmi.ft !== _static_time
        return (sim.mat_class, false)
    end

    mesh = sim.mesh
    lv = mesh.layout_version[]
    if dmi.pair_Dx !== nothing && dmi.stencil_layout == lv
        cls, R, R_ok = _fresh_mat_class(sim)
        return (cls, R_ok)
    end

    cls, R, R_ok = _fresh_mat_class(sim)
    R_ok || return (cls, false)

    vx = _class_uniform_values(sim, dmi.Dx)
    if vx === nothing
        dmi.pair_Dx = dmi.pair_Dy = dmi.pair_Dz = nothing
        dmi.stencil_layout = -1
        return (cls, false)
    end
    vy = _class_uniform_values(sim, dmi.Dy)
    if vy === nothing
        dmi.pair_Dx = dmi.pair_Dy = dmi.pair_Dz = nothing
        dmi.stencil_layout = -1
        return (cls, false)
    end
    vz = _class_uniform_values(sim, dmi.Dz)
    if vz === nothing
        dmi.pair_Dx = dmi.pair_Dy = dmi.pair_Dz = nothing
        dmi.stencil_layout = -1
        return (cls, false)
    end

    pair_x = kernel_array(_pair_table(vx))
    pair_y = kernel_array(_pair_table(vy))
    pair_z = kernel_array(_pair_table(vz))
    dmi.pair_Dx = pair_x
    dmi.pair_Dy = pair_y
    dmi.pair_Dz = pair_z

    # interfacial DMI needs a per-class D to guard the case where the current
    # cell's D_I is zero but the neighbour's D_n is not (bulk DMI inlines this
    # with `if D_I != 0`; the interfacial partition kernel can't branch on that
    # per bond, so it uses Dcls[cI+1] as the effective D_I).
    if dmi.type === :interfacial
        Dcls = zeros(T, R + 1)
        @inbounds for c in 1:R
            Dcls[c+1] = vx[c+1]   # Dx == Dy == Dz for interfacial (uniform D)
        end
        dmi.Dcls = kernel_array(Dcls)
    else
        empty!(dmi.Dcls)
    end

    dmi.stencil_layout = lv
    return (cls, true)
end
