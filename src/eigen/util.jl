# ---------------------------------------------------------------------------
# Symbolic linearisation number types (flat representation).
#
# Design rationale — replacing the old nested `AddExpr(a::AbstractFloat,
# b::AbstractFloat)` binary tree with a **single flat `Dict{Int,Float64}` of
# (ε_id → coefficient) pairs plus a scalar constant field:
#
#   value = c + Σ_{k} coefs[k] · ε_k
#
# Benefits (P0#2 + P1#5 + P2#6):
#   • O(depth) nested-tree recursion and the associated stack overflow for
#     deep expressions (N≳1000) are eliminated entirely; traversal is bounded
#     by the number of NON-ZERO ε-coefficients which is typically the 2-D
#     tangent stencil size (O(1) per cell for local-only interactions, O(N)
#     only when long-range terms were symbolically folded in — which they
#     aren't on Route Y).
#   • `collect_terms` used to allocate a fresh Dict per element call, i.e.
#     9N Dict allocations per build_matrix call; now `collect_terms` is
#     essentially a type dispatch + `copy` of the pre-existing Dict, and
#     `simplify` is the identity for already-flat types.
#   • Operator overloads perform *aggressive pruning*:
#       ε_a * ε_b = 0                   (Epsilon*Epsilon rule)
#       Real * (c + Σ…) = Real*c + Σ…   (constant-float, no AddExpr indirection)
#       0 * anything = 0                (early return, no Dict allocation)
#   • `promote_rule` is corrected: the old signatures used `::Type{Real}`
#     which is syntactically invalid because `Real` is UnionAll; promote
#     rules are only defined for CONCRETE or ABSTRACT types, not UnionAll
#     type selectors.  We define the correct signatures and keep the result
#     type `FlatTerm` for any combination that includes ε-components.
# ---------------------------------------------------------------------------

struct Epsilon <: AbstractFloat
    id::Int
    value::Float64
end

"""
    FlatTerm <: AbstractFloat

Canonical flat representation of a first-order symbolic perturbation:

    FlatTerm(c, coefs)  ≡  c + Σ_k coefs[k] · ε_k

`c` is the Float64 constant (ε⁰) term and `coefs::Dict{Int,Float64}` maps
ε-tag index → coefficient.  The Dict is kept *normalised*: entries with
abs(value) ≤ eps(Float64) are dropped, the zero-coefficient term is stored
only in `c` (never as coefs[0]), so two FlatTerms that compare equal always
have identical Dicts.
"""
struct FlatTerm <: AbstractFloat
    c::Float64
    coefs::Dict{Int,Float64}
end

# ---------------------------------------------------------------------------
# Constructors / normalisation helper
# ---------------------------------------------------------------------------

function _normalise_dict(d::Dict{Int,Float64})
    for (k, v) in d
        abs(v) <= eps(Float64) && delete!(d, k)
    end
    return d
end

function FlatTerm(c::Real, coefs::Dict{Int,<:Real})
    cc = Float64(c)
    dd = convert(Dict{Int,Float64}, copy(coefs))
    _normalise_dict(dd)
    # If only a constant survives, return a plain Float64 to keep types small.
    isempty(dd) && return cc
    return FlatTerm(cc, dd)
end

FlatTerm(x::Epsilon) = FlatTerm(0.0, Dict{Int,Float64}(x.id => x.value))
FlatTerm(x::Real)    = iszero(x) ? 0.0 : Float64(x)   # plain constant, no wrapper

# ---------------------------------------------------------------------------
# promote_rule — CORRECT signatures (no `::Type{Real}` which is UnionAll
# and invalid as a promote_rule lhs/rhs).  When mixing numeric types with
# Epsilon or FlatTerm the result is always `FlatTerm` (since any ε survives).
# ---------------------------------------------------------------------------
for (A, B, R) in [
    (:Epsilon,  :Epsilon,  :FlatTerm),
    (:FlatTerm, :Epsilon,  :FlatTerm),
    (:Epsilon,  :FlatTerm, :FlatTerm),
    (:FlatTerm, :FlatTerm, :FlatTerm),
]
    @eval Base.promote_rule(::Type{$A}, ::Type{$B}) = $R
end
for T in (:Float16, :Float32, :Float64, :Bool, :Int8, :Int16, :Int32, :Int64, :Int128,
          :UInt8, :UInt16, :UInt32, :UInt64, :UInt128, :BigInt, :BigFloat, :Rational)
    @eval Base.promote_rule(::Type{Epsilon},  ::Type{$T}) = FlatTerm
    @eval Base.promote_rule(::Type{$T},       ::Type{Epsilon})  = FlatTerm
    @eval Base.promote_rule(::Type{FlatTerm}, ::Type{$T}) = FlatTerm
    @eval Base.promote_rule(::Type{$T},       ::Type{FlatTerm}) = FlatTerm
end
# Fallbacks for any <:Integer / <:AbstractFloat / <:Rational subtypes:
Base.promote_rule(::Type{Epsilon},  ::Type{T}) where {T<:Integer}       = FlatTerm
Base.promote_rule(::Type{T},        ::Type{Epsilon})  where {T<:Integer} = FlatTerm
Base.promote_rule(::Type{FlatTerm}, ::Type{T}) where {T<:Real}          = FlatTerm
Base.promote_rule(::Type{T},        ::Type{FlatTerm}) where {T<:Real}   = FlatTerm

# ---------------------------------------------------------------------------
# Conversion / construct → promote rules above will handle most cases; we
# still need `convert` fallbacks.
# ---------------------------------------------------------------------------
Base.convert(::Type{FlatTerm}, x::Epsilon)  = FlatTerm(x)
Base.convert(::Type{FlatTerm}, x::Real)     = FlatTerm(x)
Base.convert(::Type{FlatTerm}, x::FlatTerm) = x

# AbstractFloat requires: Epsilon <: AbstractFloat → Float64(Epsilon(...)) not defined.
# Epsilon / FlatTerm pretend to be AbstractFloat but Julia's Base will try
# `convert(T, x)` when dispatching; the following stubs satisfy it:
(::Type{T})(x::Epsilon)  where {T<:AbstractFloat} = convert(T, x.value)  # used by Base kernels that do T(0) or convert from Real
(::Type{T})(x::FlatTerm) where {T<:AbstractFloat} = convert(T, x.c)       # best-effort: drop ε tags

# ---------------------------------------------------------------------------
# Multiplication — ε_a * ε_b = 0; Real scales every coefficient
# ---------------------------------------------------------------------------
Base.:*(a::Epsilon, b::Epsilon) = 0.0   # rule: ε·ε = 0 (first-order)

function Base.:*(r::Real, e::Epsilon)
    iszero(r) && return 0.0
    rv = Float64(r) * e.value
    iszero(rv) && return 0.0
    return FlatTerm(0.0, Dict{Int,Float64}(e.id => rv))
end
Base.:*(e::Epsilon, r::Real) = r * e

function Base.:*(r::Real, f::FlatTerm)
    iszero(r) && return 0.0
    rv = Float64(r)
    new_c = rv * f.c
    new_coefs = Dict{Int,Float64}()
    sizehint!(new_coefs, length(f.coefs))
    @inbounds for (k, v) in f.coefs
        nv = rv * v
        iszero(nv) && continue
        new_coefs[k] = nv
    end
    return FlatTerm(new_c, new_coefs)   # will normalise
end
Base.:*(f::FlatTerm, r::Real) = r * f

function Base.:*(e::Epsilon, f::FlatTerm)
    # Epsilon carries a single ε term; multiplying into FlatTerm must apply
    # ε*ε = 0 to every ε term inside f → only the constant survives scaling.
    return e * f.c
end
Base.:*(f::FlatTerm, e::Epsilon) = e * f

function Base.:*(a::FlatTerm, b::FlatTerm)
    # FOIL, but any product containing two ε indices vanishes:
    #   (a0 + Σ a_k ε_k) * (b0 + Σ b_j ε_j)
    # = a0*b0 + a0*(Σ b_j ε_j) + (Σ a_k ε_k)*b0   [cross terms with two ε vanish]
    # = a0*b0 + Σ (a0*b_j + b0*a_k) ε_j [with k/j merged]
    c_new = a.c * b.c
    out = Dict{Int,Float64}()
    if !iszero(a.c)
        ac = a.c
        for (k, v) in b.coefs
            nv = ac * v
            iszero(nv) && continue
            out[k] = get(out, k, 0.0) + nv
        end
    end
    if !iszero(b.c)
        bc = b.c
        for (k, v) in a.coefs
            nv = bc * v
            iszero(nv) && continue
            out[k] = get(out, k, 0.0) + nv
        end
    end
    return FlatTerm(c_new, out)
end

# ---------------------------------------------------------------------------
# Addition / subtraction — pointwise Dict merge with coefficient summation
# ---------------------------------------------------------------------------
function Base.:+(a::Epsilon, b::Epsilon)
    if a.id == b.id
        s = a.value + b.value
        iszero(s) && return 0.0
        return FlatTerm(0.0, Dict{Int,Float64}(a.id => s))
    end
    return FlatTerm(0.0, Dict{Int,Float64}(a.id => a.value, b.id => b.value))
end

function Base.:+(r::Real, e::Epsilon)
    iszero(r) && return FlatTerm(e)
    return FlatTerm(Float64(r), Dict{Int,Float64}(e.id => e.value))
end
Base.:+(e::Epsilon, r::Real) = r + e

function Base.:+(r::Real, f::FlatTerm)
    iszero(r) && return f
    return FlatTerm(f.c + Float64(r), f.coefs)
end
Base.:+(f::FlatTerm, r::Real) = r + f

function Base.:+(e::Epsilon, f::FlatTerm)
    d = copy(f.coefs)
    d[e.id] = get(d, e.id, 0.0) + e.value
    return FlatTerm(f.c, d)   # normalises
end
Base.:+(f::FlatTerm, e::Epsilon) = e + f

function Base.:+(a::FlatTerm, b::FlatTerm)
    # Merge b's coefs into a copy of a.coefs
    if length(a.coefs) >= length(b.coefs)
        d = copy(a.coefs)
        for (k, v) in b.coefs
            d[k] = get(d, k, 0.0) + v
        end
    else
        d = copy(b.coefs)
        for (k, v) in a.coefs
            d[k] = get(d, k, 0.0) + v
        end
    end
    return FlatTerm(a.c + b.c, d)   # normalises
end

# Unary minus
Base.:-(e::Epsilon)  = Epsilon(e.id, -e.value)
Base.:-(f::FlatTerm) = FlatTerm(-f.c, Dict{Int,Float64}(k => -v for (k, v) in f.coefs))

# Subtraction: reuse addition of negation
Base.:-(a::Epsilon, b::Epsilon)  = a + (-b)
Base.:-(r::Real,    e::Epsilon)  = r + (-e)
Base.:-(e::Epsilon, r::Real)     = e + (-r)
Base.:-(r::Real,    f::FlatTerm) = r + (-f)
Base.:-(f::FlatTerm, r::Real)    = f + (-r)
Base.:-(e::Epsilon, f::FlatTerm) = e + (-f)
Base.:-(f::FlatTerm, e::Epsilon) = f + (-e)
Base.:-(a::FlatTerm, b::FlatTerm) = a + (-b)

# ---------------------------------------------------------------------------
# Pretty printing — mirrors the old format (e.g. "0.5 + 1.2ε_3 - 0.7ε_1"),
# so the user-facing experience doesn't change.
# ---------------------------------------------------------------------------
function Base.show(io::IO, x::AbstractFloat)
    if x isa Epsilon
        if x.value >= 0
            print(io, x.value, "ε_", x.id)
        else
            print(io, "(", x.value, "ε_", x.id, ")")
        end
    elseif x isa FlatTerm
        print(io, "(")
        first = true
        if !iszero(x.c) || isempty(x.coefs)
            print(io, x.c)
            first = false
        end
        for k in sort(collect(keys(x.coefs)))
            v = x.coefs[k]
            iszero(v) && continue
            if first
                if v >= 0
                    print(io, v, "ε_", k)
                else
                    print(io, "(", v, "ε_", k, ")")
                end
                first = false
            else
                if v >= 0
                    print(io, " + ", v, "ε_", k)
                else
                    print(io, " - ", -v, "ε_", k)
                end
            end
        end
        print(io, ")")
    else
        print(io, x)
    end
end

# ---------------------------------------------------------------------------
# collect_terms — returns a Dict{Int,Float64}: id=0 is the constant, id≥1 are
# the ε coefficients.  Keeps the public contract used by eigen.jl's B-matrix
# extraction stage.  For FlatTerm this is O(nonzero) + small Dict allocation.
# For Epsilon it's O(1).  For plain Real it's two-element Dict.
# ---------------------------------------------------------------------------
function collect_terms(x::AbstractFloat)
    if x isa FlatTerm
        d = Dict{Int,Float64}()
        sizehint!(d, 1 + length(x.coefs))
        iszero(x.c) || (d[0] = x.c)
        for (k, v) in x.coefs
            iszero(v) && continue
            d[k] = v
        end
        return d
    elseif x isa Epsilon
        return Dict{Int,Float64}(x.id => x.value)
    else  # plain Real → constant term
        f = Float64(x)
        return iszero(f) ? Dict{Int,Float64}() : Dict{Int,Float64}(0 => f)
    end
end

# ---------------------------------------------------------------------------
# simplify — already-flat types come back canonical (constants promoted to
# Float64 when the Dict is empty; zero-coef entries dropped during
# construction).  For legacy callers that still pass a nested-AddExpr tree
# (shouldn't happen any more but we keep it safe), `collect_terms` flattens
# it and we re-wrap into the canonical form.
# ---------------------------------------------------------------------------
function simplify(x::AbstractFloat)
    if x isa FlatTerm
        isempty(x.coefs) && return x.c
        return x
    elseif x isa Epsilon
        return x
    elseif x isa Real   # plain Float/Int — pass through
        return Float64(x)
    else
        # Unknown symbolic subtype (shouldn't be hit once AddExpr is gone):
        # flatten via collect_terms.
        t = collect_terms(x)
        c = get(t, 0, 0.0)
        delete!(t, 0)
        return FlatTerm(c, t)
    end
end

# ---------------------------------------------------------------------------
# Compatibility alias — the old `AddExpr` name is no longer used but
# `eigen.jl` refers to it in `SymT = Union{…}`.  We keep it pointing at
# FlatTerm so callers that still mention `AddExpr` (type annotations, Union
# membership) silently resolve to the new representation.  Callers that
# manually construct `AddExpr(a,b)` nodes are intentionally broken — that
# tree-building pattern is exactly what we're eliminating.
# ---------------------------------------------------------------------------
const AddExpr = FlatTerm

# ---------------------------------------------------------------------------
# Rotation matrices — numerical stability + precomputation hooks.
#
# P2 #4 / #8 — the original rotation_matrix(mx,my,mz) used atan / acos on the
# direction cosines.  This has two issues:
#   1. Numerical conditioning when |m_z| ≈ 1 (tan φ = y/x → 0/0 singularity
#      when m is near ±ẑ; θ near 0 or π forces sin θ cancellation).
#   2. The inverse `rotation_matrix_inverse` calls `rotation_matrix` again
#      and then transposes — duplicating the trig work.  For N large the
#      duplicates add up; eigen.jl calls R *and* R_inv per cell *per* stage.
#
# We replace both with a single numerically-stable orthogonal construction:
# given unit vector m, build a right-handed orthonormal basis (q1, q2, m)
# using a Householder-type reflection step that avoids the sin θ = 0
# singularity and uses only sqrt, div, mul (no trig functions at all).
# The inverse is still the transpose and is returned from the same internal
# helper when both are requested.  For performance the public API remains
# `rotation_matrix` / `rotation_matrix_inverse` (returning fresh Matrix)
# because most callers expect that; we also expose in-place helpers that the
# eigen path can use to avoid 3×3 heap allocations per cell.
# ---------------------------------------------------------------------------

function _ortho_basis_q1q2!(Q, mx, my, mz)
    # Given unit (mx,my,mz) fill Q[1:3,1] and Q[1:3,2] with orthonormal
    # vectors orthogonal to (mx,my,mz), forming right-handed basis (q1,q2,m̂).
    # Uses: choose the component of m̂ with smallest absolute value for the
    # Householder seed — avoids catastrophic cancellation — see
    # "Building an Orthonormal Basis from a Unit Vector", Duff et al. 2017.
    abs_x, abs_y, abs_z = abs(mx), abs(my), abs(mz)
    if abs_x <= abs_y && abs_x <= abs_z
        # smallest component is x → (0, -mz, my)
        a = my*my + mz*mz   # = 1 - mx^2 (guaranteed non-zero if mx smallest and |m|=1)
        s = sqrt(a)
        qx1, qy1, qz1 = 0.0, -mz/s, my/s
        # q2 = m × q1  (right-handed)
        qx2 = my*qz1 - mz*qy1
        qy2 = mz*qx1 - mx*qz1
        qz2 = mx*qy1 - my*qx1
    elseif abs_y <= abs_z
        # smallest component is y → (mz, 0, -mx)
        a = mx*mx + mz*mz
        s = sqrt(a)
        qx1, qy1, qz1 = mz/s, 0.0, -mx/s
        qx2 = my*qz1 - mz*qy1
        qy2 = mz*qx1 - mx*qz1
        qz2 = mx*qy1 - my*qx1
    else
        # smallest component is z → (-my, mx, 0)
        a = mx*mx + my*my
        s = sqrt(a)
        qx1, qy1, qz1 = -my/s, mx/s, 0.0
        qx2 = my*qz1 - mz*qy1
        qy2 = mz*qx1 - mx*qz1
        qz2 = mx*qy1 - my*qx1
    end
    Q[1,1] = qx1;  Q[2,1] = qy1;  Q[3,1] = qz1
    Q[1,2] = qx2;  Q[2,2] = qy2;  Q[3,2] = qz2
    Q[1,3] = mx;   Q[2,3] = my;   Q[3,3] = mz
    return Q
end

raw"""
R z = m0   (R = orthonormal, right-handed)

Columns of R are (tangent-1, tangent-2, m̂0).  Uses a singularity-free
Householder-type basis construction (no atan / acos, no conditioning issues
when m̂0 is near ±ẑ).
"""
function rotation_matrix(mx, my, mz)
    n = sqrt(mx*mx + my*my + mz*mz)
    if n <= eps(Float64)
        # Degenerate zero-length magnetisation: identity.
        A = zeros(3, 3)
        A[1,1] = 1;  A[2,2] = 1;  A[3,3] = 1
        return A
    end
    # Normalise once — downstream logic assumes unit length.
    Q = Matrix{Float64}(undef, 3, 3)
    _ortho_basis_q1q2!(Q, mx/n, my/n, mz/n)
    return Q
end

"""
R^{-1} = R^{T} because R is orthogonal.  Same stable basis construction
as `rotation_matrix`, but returns the transpose directly (so we don't
allocate R just to transpose it).
"""
function rotation_matrix_inverse(mx, my, mz)
    n = sqrt(mx*mx + my*my + mz*mz)
    if n <= eps(Float64)
        A = zeros(3, 3)
        A[1,1] = 1;  A[2,2] = 1;  A[3,3] = 1
        return A
    end
    Q = Matrix{Float64}(undef, 3, 3)
    _ortho_basis_q1q2!(Q, mx/n, my/n, mz/n)
    # In-place transpose would be nice, but Q is only 3×3; permute explicitly.
    return @inbounds [
        Q[1,1]  Q[2,1]  Q[3,1];
        Q[1,2]  Q[2,2]  Q[3,2];
        Q[1,3]  Q[2,3]  Q[3,3];
    ]
end

# ---------------------------------------------------------------------------
# 3-D cross product components — already optimal; kept for API stability.
# ---------------------------------------------------------------------------
@inline function cross_x(x1, x2, x3, y1, y2, y3)
    return -x3 * y2 + x2 * y3
end
@inline function cross_y(x1, x2, x3, y1, y2, y3)
    return x3 * y1 - x1 * y3
end
@inline function cross_z(x1, x2, x3, y1, y2, y3)
    return -x2 * y1 + x1 * y2
end
