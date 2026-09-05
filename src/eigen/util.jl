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

    # Public constructor: copy + exact-zero prune + collapse to plain Float64
    # when no ε survives — semantics unchanged (moved INTO the struct: outer
    # constructors cannot call `new`, and the previous outer-then-inner chain
    # copied every Dict twice; SFLAT_PLAN.md §2.1).
    function FlatTerm(c::Real, coefs::Dict{Int,<:Real})
        cc = Float64(c)
        dd = convert(Dict{Int,Float64}, copy(coefs))
        _normalise_dict(dd)
        isempty(dd) && return cc
        return new(cc, dd)
    end

    # Raw constructor: trust (c, d) — no copy, no prune, no collapse.  For the
    # arithmetic ops that just built a fresh Dict (kills the second copy).
    FlatTerm(c::Float64, d::Dict{Int,Float64}, ::Val{:raw}) = new(c, d)
end

# ---------------------------------------------------------------------------
# Constructors / normalisation helper
# ---------------------------------------------------------------------------

function _normalise_dict(d::Dict{Int,Float64})
    # Exact-zero pruning only (SFLAT_PLAN.md §2.1): the old abs(v) ≤ eps(Float64)
    # scan dropped ~1e-16 dust, but it ran on EVERY op construction. Downstream
    # extraction thresholds at abs(v) > eps(Float64) anyway, so keeping dust in
    # the Dict changes no downstream result while saving one full scan per op.
    for (k, v) in d
        iszero(v) && delete!(d, k)
    end
    return d
end

FlatTerm(x::Epsilon) = FlatTerm(0.0, Dict{Int,Float64}(x.id => x.value))
FlatTerm(x::Real)    = iszero(x) ? 0.0 : Float64(x)   # plain constant, no wrapper

# Helper for the arithmetic ops: they build a FRESH Dict and previously handed
# it to the public constructor, which copied it a SECOND time — one Dict copy
# per op wasted.  `_flatterm` skips the copy/prune (raw inner constructor) and
# keeps the empty-Dict → plain-Float64 collapse.  Ops that reuse an EXISTING
# Dict (e.g. `Real + FlatTerm` reusing f.coefs) must keep the copying
# constructor for alias safety.
@inline _flatterm(c::Real, d::Dict{Int,Float64}) =
    isempty(d) ? Float64(c) : FlatTerm(Float64(c), d, Val(:raw))

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
# Epsilon / FlatTerm pretend to be AbstractFloat, but we REFUSE to silently
# drop ε via `(::Type{T})(::Epsilon/FlatTerm) where T<:AbstractFloat`. The
# previous stubs (returning `x.value` / `x.c`) were the enabler of the
# KernelAbstractions CPU-backend silent-truncation bug: when a kernel
# accumulator was lowered to a `Float64` stack slot, any `+=` of an
# Epsilon/FlatTerm value went through `Float64(::FlatTerm)` and quietly
# discarded every ε coefficient.  By throwing here we make every kernel
# swallow point LOUD (MethodError), so the only way to extract the constant
# part is the explicit `cpart(x)` accessor below — code that genuinely
# wants to drop ε must ask for it by name.
function (::Type{T})(x::Union{Epsilon,FlatTerm}) where {T<:AbstractFloat}
    error("Refusing to convert $(typeof(x)) to $T — this would silently drop ε. " *
          "Use cpart(x) for the constant term, or zero(T)/one(T) for type-correct init.")
end

"""
    cpart(x) -> Float64

Return the Float64 constant (ε⁰) part of a symbolic number.  Use this when
you genuinely want to discard the ε coefficients (e.g. extracting a
Float64 baseline from a symbolic field buffer).  For symbolic numbers
stored in `Vector{AbstractFloat}` this is the safe way to materialise a
`Vector{Float64}` without triggering the loud `convert` stub above.
"""
cpart(x::Real)       = Float64(x)
cpart(x::Epsilon)    = 0.0
cpart(x::FlatTerm)   = x.c
cpart(xs::AbstractArray) = map(cpart, xs)
export cpart

# ---------------------------------------------------------------------------
# Multiplication — ε_a * ε_b = 0; Real scales every coefficient
# ---------------------------------------------------------------------------
Base.:*(a::Epsilon, b::Epsilon) = 0.0   # rule: ε·ε = 0 (first-order)

function Base.:*(r::Real, e::Epsilon)
    iszero(r) && return 0.0
    rv = Float64(r) * e.value
    iszero(rv) && return 0.0
    return _flatterm(0.0, Dict{Int,Float64}(e.id => rv))
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
    return _flatterm(new_c, new_coefs)
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
    return _flatterm(c_new, out)
end

# ---------------------------------------------------------------------------
# zero / one — CRITICAL overrides.  Without these, `zero(FlatTerm)` falls
# back to `convert(FlatTerm, 0)` → FlatTerm(::Int<:Real) at line 76 returns
# a bare `0.0 :: Float64`, which strips ALL ε-coefficients from fx/fy/fz in
# the kernel immediately.  Same pitfall for Epsilon.
# ---------------------------------------------------------------------------
Base.zero(::Type{Epsilon})  = FlatTerm(0.0, Dict{Int,Float64}())
Base.zero(::Epsilon)        = zero(Epsilon)
Base.zero(::Type{FlatTerm}) = FlatTerm(0.0, Dict{Int,Float64}())
Base.zero(f::FlatTerm)      = FlatTerm(0.0, Dict{Int,Float64}())   # don't reuse f.coefs

Base.one(::Type{Epsilon})   = FlatTerm(1.0, Dict{Int,Float64}())
Base.one(::Epsilon)         = one(Epsilon)
Base.one(::Type{FlatTerm})  = FlatTerm(1.0, Dict{Int,Float64}())
Base.one(f::FlatTerm)       = FlatTerm(1.0, Dict{Int,Float64}())

# ---------------------------------------------------------------------------
# Symbolic-aware cross product — inline the 3 explicit component formulas
# so NO promote/convert EVER touches the mixed Float64 × SymType arguments.
# This avoids the ambiguity chain caused by util.jl's generic fallback
# promoting 6 args to a common type before compute, which strips ε tags.
# Note: cross_x/y/z are exported from MicroMagnetic, so we add method
# overloads that dispatch whenever ANY of the 6 args is Epsilon/FlatTerm.
# ---------------------------------------------------------------------------
const SymType = Union{Epsilon, FlatTerm}

# ---------------------------------------------------------------------------
# Division / inversion.  Julia's generic `/(x, y)` falls back to
# `promote(x, y)` first, which throws `promotion of types FlatTerm & Float64
# failed to change any arguments` because we deliberately kept FlatTerm ≠
# Float64 for symbolic clarity.  Add explicit dispatches so:
#
#   FlatTerm / Real  →  FlatTerm * inv(Real)  (trivial, no ε terms)
#   Real / FlatTerm  →  Real * inv(FlatTerm)
#   inv(FlatTerm)    →  allowed ONLY when the symbolic part is empty, i.e.
#                       the FlatTerm represents a bare constant.  Otherwise
#                       we cannot divide 1/(c + Σcoef·ε) at FIRST ORDER.
# ---------------------------------------------------------------------------
Base.inv(e::Epsilon)  = error("Cannot invert Epsilon perturbation: 1/ε is not first-order.")

function Base.inv(f::FlatTerm)
    if isempty(f.coefs)
        # Pure constant: just 1/c.
        invc = inv(f.c)
        iszero(invc) && return 0.0
        return invc   # Plain Float64 since there's no ε part.
    end
    # 1/(c + Σ) at first order:  1/c · (1 − Σ/c)  →  return as FlatTerm.
    # If c == 0 here (singular but has coefs), it's 1/ε-class: error out.
    iszero(f.c) && error("Cannot invert FlatTerm with zero constant and non-zero ε part — not first-order representable.")
    invc = inv(f.c)
    new_coefs = Dict{Int,Float64}()
    sizehint!(new_coefs, length(f.coefs))
    coef = -invc * invc
    for (k, v) in f.coefs
        nv = coef * v
        iszero(nv) && continue
        new_coefs[k] = nv
    end
    return _flatterm(invc, new_coefs)
end

# Left-division too:  `Real \ SymType`  =  inv(Real) * SymType.
Base.:/(a::SymType,    b::Real)     = a * inv(b)
Base.:\(a::Real,       b::SymType)  = inv(a) * b
Base.:/(a::Epsilon,    b::Epsilon)  = error("Cannot divide Epsilon/Epsilon at first order (ε/ε is O(1) indeterminate).")
Base.:/(a::Epsilon,    b::FlatTerm) = a * inv(b)
Base.:/(a::FlatTerm,   b::Epsilon)  = error("Cannot divide FlatTerm / Epsilon at first order (divide by ε is singular).")
Base.:/(a::FlatTerm,   b::FlatTerm) = a * inv(b)
Base.:/(a::Real,       b::SymType)  = a * inv(b)

# ---------------------------------------------------------------------------
# Addition / subtraction — pointwise Dict merge with coefficient summation
# ---------------------------------------------------------------------------
function Base.:+(a::Epsilon, b::Epsilon)
    if a.id == b.id
        s = a.value + b.value
        iszero(s) && return 0.0
        return _flatterm(0.0, Dict{Int,Float64}(a.id => s))
    end
    return _flatterm(0.0, Dict{Int,Float64}(a.id => a.value, b.id => b.value))
end

function Base.:+(r::Real, e::Epsilon)
    iszero(r) && return FlatTerm(e)
    return _flatterm(Float64(r), Dict{Int,Float64}(e.id => e.value))
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
    return _flatterm(f.c, d)   # normalises
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
    return _flatterm(a.c + b.c, d)   # normalises
end

# Unary minus
Base.:-(e::Epsilon)  = Epsilon(e.id, -e.value)
Base.:-(f::FlatTerm) = _flatterm(-f.c, Dict{Int,Float64}(k => -v for (k, v) in f.coefs))

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

@doc raw"""
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
# 3-D cross product components  cross_x / cross_y / cross_z.
#
# DEFINITION SOURCE:  src/util.jl  (NOT here!)
#   As of 2026-08-29 util.jl provides BOTH:
#     (a) a fully generic, untyped  cross_x(x1,x2,x3,y1,y2,y3)  fallback that
#         accepts mixed-type argument tuples like (Float64,Float64,Float64,
#         AddExpr,AddExpr,AddExpr) which the symbolic build_matrix linearisation
#         relies on for DMI / anisotropy / vector cross kernels;
#     (b) the  cross_x(::T,...) where {T<:Number}  homogeneous fast-path used by
#         the production Float64 simulation kernels (fully inlined).
#
# DO NOT RE-DEFINE cross_x / cross_y / cross_z here.  Julia treats any method
# with the same signature in a later include() as an *overwrite*, which is
# FORBIDDEN during module precompilation ("ERROR: Method overwriting is not
# permitted during Module precompilation").  Precompilation used to silently
# load a half-recompiled MicroMagnetic.jl image which in turn caused the
# bulkdmi_kernel! DMI linearisation path to return all zeros even for real-
# valued inputs; deleting this duplicate block fixed both the precompilation
# blocker and the downstream DMI silent-zero bug.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# SFlat{CAP} — fixed-capacity isbits symbolic dual numbers (SFLAT_PLAN.md §2.2)
#
#   value = c + Σ_{j≤n} vs[j] · ε_{ids[j]}
#
# Same first-order semantics as FlatTerm (ε·ε = 0, same-id merge, constant
# collapse), but the ε coefficients live in inline NTuples instead of a heap
# Dict: isbits values, zero allocation per op, concrete-typed buffers (no
# boxing).  Prototype measurements (2026-09-05, /tmp/fd_bench/sflat_proto.jl):
# single-tag 82× vs FlatTerm (parity with ForwardDiff.Dual{Tag,Float64,1}),
# dense multi-tag 10×, results bitwise identical in both regimes.
#
# CAP is a compile-time bound on the number of DISTINCT tags a scalar may
# carry; exceeding it throws loudly (choose a larger CAP via
# set_precision(SFlat{CAP}), or fall back to the unbounded
# set_precision(AbstractFloat) mode).  3-D six-neighbour stencils stay ≤14
# tags, so CAP = 16 covers exchange/DMI/anisotropy plus the LLG cross terms.
#
# `<: AbstractFloat` mirrors the Epsilon trick: existing kernel signatures
# `where {T<:AbstractFloat}` accept SFlat unchanged.  Kernel bodies are
# number-type clean — verified bitwise against the Float64 path.
# ---------------------------------------------------------------------------
struct SFlat{CAP} <: AbstractFloat
    c::Float64
    n::Int                            # active tag count (≤ CAP)
    ids::NTuple{CAP,Int32}
    vs::NTuple{CAP,Float64}
end

@inline _sflat_zids(::Val{CAP}) where CAP = ntuple(_ -> Int32(0), Val(CAP))
@inline _sflat_zvs(::Val{CAP})  where CAP = ntuple(_ -> 0.0, Val(CAP))

Base.zero(::Type{SFlat{CAP}}) where CAP =
    SFlat{CAP}(0.0, 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))
Base.one(::Type{SFlat{CAP}}) where CAP =
    SFlat{CAP}(1.0, 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))
SFlat{CAP}(x::Real) where CAP =
    SFlat{CAP}(Float64(x), 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))

# Mixing symbolic modes is a bug (tags would be silently dropped) — loud,
# mirroring the Epsilon/FlatTerm convert stub above.
SFlat{CAP}(x::Union{Epsilon,FlatTerm}) where CAP =
    error("Refusing to convert $(typeof(x)) to SFlat{$CAP} — this would silently drop ε. " *
          "Use one symbolic mode per sim (set_precision(AbstractFloat) or set_precision(SFlat{CAP})).")

Base.convert(::Type{SFlat{CAP}}, x::Real) where CAP = SFlat{CAP}(x)
Base.promote_rule(::Type{SFlat{CAP}}, ::Type{SFlat{CAP}}) where CAP = SFlat{CAP}
Base.promote_rule(::Type{SFlat{CAP}}, ::Type{T}) where {CAP, T<:Real} = SFlat{CAP}
Base.promote_rule(::Type{T}, ::Type{SFlat{CAP}}) where {CAP, T<:Real} = SFlat{CAP}
# Cross-mode / cross-CAP promotion falls back to the unbounded FlatTerm.
Base.convert(::Type{FlatTerm}, x::SFlat) = FlatTerm(x.c, collect_terms(x))
Base.promote_rule(::Type{SFlat{CAP}}, ::Type{Epsilon})  where CAP = FlatTerm
Base.promote_rule(::Type{Epsilon},  ::Type{SFlat{CAP}}) where CAP = FlatTerm
Base.promote_rule(::Type{SFlat{CAP}}, ::Type{FlatTerm}) where CAP = FlatTerm
Base.promote_rule(::Type{FlatTerm}, ::Type{SFlat{CAP}}) where CAP = FlatTerm
Base.promote_rule(::Type{SFlat{A}}, ::Type{SFlat{B}}) where {A, B} = FlatTerm

# Constant-only SFlat values (the kind stored in parameter buffers) convert
# losslessly to plain floats — this is what lets the Float64 demag path read
# SFlat-mode parameter buffers.  Tagged values refuse loudly, same protection
# as the Epsilon/FlatTerm stub above (KA CPU silent-truncation guard).
for Tf in (:Float16, :Float32, :Float64)
    @eval function Base.convert(::Type{$Tf}, x::SFlat)
        iszero(x.n) || error("Refusing to convert a tagged SFlat to $($Tf) — this would " *
                             "silently drop ε. Use cpart(x) / coefof(x, id).")
        return $Tf(x.c)
    end
end
Base.Float64(x::SFlat) = convert(Float64, x)

@inline function _sflat_insert(n::Int, ids::NTuple{CAP,Int32}, vs::NTuple{CAP,Float64},
                               id::Int32, v::Float64) where {CAP}
    @inbounds for t in 1:n
        ids[t] == id && return n, ids, Base.setindex(vs, vs[t] + v, t)
    end
    n + 1 <= CAP || throw(OverflowError("SFlat: tag capacity $CAP exceeded (need $(n + 1)). " *
        "Use set_precision(SFlat{$(n + 2)}) or larger, or the unbounded set_precision(AbstractFloat) mode."))
    return n + 1, Base.setindex(ids, id, n + 1), Base.setindex(vs, v, n + 1)
end

@inline function Base.:+(a::SFlat{CAP}, b::SFlat{CAP}) where {CAP}
    # Constant fast paths: SFlat-mode parameter buffers hold n = 0 values, so
    # bond arithmetic mostly hits these.
    iszero(a.n) && return SFlat{CAP}(a.c + b.c, b.n, b.ids, b.vs)
    iszero(b.n) && return SFlat{CAP}(a.c + b.c, a.n, a.ids, a.vs)
    n, ids, vs = a.n, a.ids, a.vs
    @inbounds for j in 1:b.n
        n, ids, vs = _sflat_insert(n, ids, vs, b.ids[j], b.vs[j])
    end
    return SFlat{CAP}(a.c + b.c, n, ids, vs)
end

@inline _sflat_scale(r::Float64, a::SFlat{CAP}) where CAP =
    SFlat{CAP}(r * a.c, a.n, a.ids, map(v -> r * v, a.vs))

@inline _sflat_const(::Val{CAP}, v::Float64) where CAP =
    SFlat{CAP}(v, 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))

# ---------------------------------------------------------------------------
# First-order primitive rules (ForwardDiff-style, minimal flops).  For a
# constant-only SFlat (n = 0) these reduce to the plain function applied to
# the constant — this is what makes host-side helpers work in SFlat mode
# (init_m0's normalise kernel calls sqrt on buffer values, and mesh
# coordinates are precision-typed, so user-supplied m0 functions receive
# SFlat constants).
# ---------------------------------------------------------------------------
@inline function Base.sqrt(a::SFlat{CAP}) where {CAP}
    c = sqrt(a.c)
    iszero(a.n) && return _sflat_const(Val(CAP), c)
    iszero(a.c) &&
        error("SFlat sqrt: derivative of sqrt at a zero constant is singular.")
    s = 1.0 / (2.0 * c)
    return SFlat{CAP}(c, a.n, a.ids, map(v -> s * v, a.vs))
end

@inline function Base.sin(a::SFlat{CAP}) where {CAP}
    c = sin(a.c)
    iszero(a.n) && return _sflat_const(Val(CAP), c)
    s = cos(a.c)
    return SFlat{CAP}(c, a.n, a.ids, map(v -> s * v, a.vs))
end

@inline function Base.cos(a::SFlat{CAP}) where {CAP}
    c = cos(a.c)
    iszero(a.n) && return _sflat_const(Val(CAP), c)
    s = -sin(a.c)
    return SFlat{CAP}(c, a.n, a.ids, map(v -> s * v, a.vs))
end

@inline function Base.exp(a::SFlat{CAP}) where {CAP}
    c = exp(a.c)
    iszero(a.n) && return _sflat_const(Val(CAP), c)
    return SFlat{CAP}(c, a.n, a.ids, map(v -> c * v, a.vs))
end

@inline function Base.log(a::SFlat{CAP}) where {CAP}
    c = log(a.c)          # DomainError for c ≤ 0 comes from Base naturally
    iszero(a.n) && return _sflat_const(Val(CAP), c)
    s = 1.0 / a.c
    return SFlat{CAP}(c, a.n, a.ids, map(v -> s * v, a.vs))
end

@inline function Base.:^(a::SFlat{CAP}, p::Integer) where {CAP}
    iszero(a.n) && return _sflat_const(Val(CAP), a.c^p)
    c = a.c^p
    s = p * a.c^(p - 1)
    return SFlat{CAP}(c, a.n, a.ids, map(v -> s * v, a.vs))
end

@inline function Base.:^(a::SFlat{CAP}, p::Real) where {CAP}
    iszero(a.n) && return _sflat_const(Val(CAP), a.c^p)
    a.c > 0 ||
        error("SFlat ^: non-integer power requires a positive constant part (got $(a.c)).")
    c = a.c^p
    s = p * a.c^(p - 1)
    return SFlat{CAP}(c, a.n, a.ids, map(v -> s * v, a.vs))
end

@inline Base.:^(r::Real, a::SFlat{CAP}) where CAP = iszero(a.n) ?
    _sflat_const(Val(CAP), Float64(r)^a.c) : exp(a * log(Float64(r)))

@inline function Base.:*(a::SFlat{CAP}, b::SFlat{CAP}) where {CAP}
    (iszero(a.c) && iszero(a.n)) && return zero(SFlat{CAP})
    (iszero(b.c) && iszero(b.n)) && return zero(SFlat{CAP})
    # ε·ε = 0 (first order): only constant×tag products survive — same rule as FlatTerm.
    iszero(a.n) && iszero(b.n) &&
        return SFlat{CAP}(a.c * b.c, 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))
    iszero(a.n) && return _sflat_scale(a.c, b)
    iszero(b.n) && return _sflat_scale(b.c, a)
    n, ids, vs = 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP))
    if !iszero(b.c)
        @inbounds for j in 1:a.n
            n, ids, vs = _sflat_insert(n, ids, vs, a.ids[j], a.vs[j] * b.c)
        end
    end
    if !iszero(a.c)
        @inbounds for j in 1:b.n
            n, ids, vs = _sflat_insert(n, ids, vs, b.ids[j], b.vs[j] * a.c)
        end
    end
    return SFlat{CAP}(a.c * b.c, n, ids, vs)
end

@inline Base.:*(r::Real, a::SFlat{CAP}) where CAP =
    iszero(r) ? zero(SFlat{CAP}) : _sflat_scale(Float64(r), a)
@inline Base.:*(a::SFlat{CAP}, r::Real) where CAP = r * a
@inline Base.:+(r::Real, a::SFlat{CAP}) where CAP = SFlat{CAP}(r + a.c, a.n, a.ids, a.vs)
@inline Base.:+(a::SFlat{CAP}, r::Real) where CAP = r + a
@inline Base.:-(a::SFlat{CAP}) where CAP =
    SFlat{CAP}(-a.c, a.n, a.ids, map(v -> -v, a.vs))
@inline Base.:-(a::SFlat{CAP}, b::SFlat{CAP}) where CAP = a + (-b)
@inline Base.:-(r::Real, a::SFlat{CAP}) where CAP = r + (-a)
@inline Base.:-(a::SFlat{CAP}, r::Real) where CAP = a + (-r)

# First-order inversion, mirroring inv(::FlatTerm): 1/(c + Σ) = 1/c − Σ/c².
@inline function Base.inv(a::SFlat{CAP}) where {CAP}
    iszero(a.n) && return SFlat{CAP}(1.0 / a.c, 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))
    iszero(a.c) && error("Cannot invert an SFlat with zero constant and non-zero ε part — " *
                         "not first-order representable.")
    invc = 1.0 / a.c
    n, ids, vs = 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP))
    @inbounds for j in 1:a.n
        n, ids, vs = _sflat_insert(n, ids, vs, a.ids[j], -invc * invc * a.vs[j])
    end
    return SFlat{CAP}(invc, n, ids, vs)
end
Base.:/(a::SFlat{CAP}, b::SFlat{CAP}) where CAP = a * inv(b)
@inline Base.:/(a::SFlat{CAP}, r::Real) where CAP = a * inv(Float64(r))

# Comparisons act on the constant part (ForwardDiff-style value semantics).
@inline Base.:(==)(a::SFlat, b::Real)   = a.c == b
@inline Base.:(==)(a::Real, b::SFlat)   = b.c == a
@inline Base.:(==)(a::SFlat, b::SFlat)  = a.c == b.c
@inline Base.:<(a::SFlat, b::Real)      = a.c < b
@inline Base.:<(a::Real, b::SFlat)      = a < b.c
@inline Base.:<(a::SFlat, b::SFlat)     = a.c < b.c
@inline Base.:<=(a::SFlat, b::SFlat)    = a.c <= b.c
@inline Base.isless(a::SFlat, b::SFlat) = isless(a.c, b.c)
@inline Base.isless(a::SFlat, b::Real)  = isless(a.c, b)
@inline Base.isless(a::Real, b::SFlat)  = isless(a, b.c)

cpart(x::SFlat) = x.c

function collect_terms(x::SFlat)
    d = Dict{Int,Float64}()
    sizehint!(d, 1 + x.n)
    iszero(x.c) || (d[0] = x.c)
    @inbounds for j in 1:x.n
        iszero(x.vs[j]) && continue
        d[Int(x.ids[j])] = x.vs[j]
    end
    return d
end

simplify(x::SFlat) = iszero(x.n) ? x.c : x

"""
    coefof(x, id) -> Float64

The ε_id coefficient of a symbolic (or plain) number — a direct field/Dict
lookup with NO per-call Dict copy, unlike `collect_terms`.  Hot-path accessor
for the eigen extraction stages (`_δH_loc_single!`, dense B assembly).
"""
coefof(x::Real, id::Int)      = 0.0
coefof(x::Epsilon, id::Int)   = x.id == id ? x.value : 0.0
coefof(x::FlatTerm, id::Int)  = get(x.coefs, id, 0.0)
function coefof(x::SFlat, id::Int)
    @inbounds for j in 1:x.n
        x.ids[j] == id && return x.vs[j]
    end
    return 0.0
end

@inline function mkterm(::Type{SFlat{CAP}}, c::Real, id::Int, v::Real) where CAP
    iszero(v) && return SFlat{CAP}(Float64(c), 0, _sflat_zids(Val(CAP)), _sflat_zvs(Val(CAP)))
    return SFlat{CAP}(Float64(c), 1,
                      ntuple(i -> i == 1 ? Int32(id) : Int32(0), Val(CAP)),
                      ntuple(i -> i == 1 ? Float64(v) : 0.0, Val(CAP)))
end
@inline function mkterm(::Type{SFlat{CAP}}, c::Real, id1::Int, v1::Real, id2::Int, v2::Real) where CAP
    n1 = !iszero(v1)
    n2 = !iszero(v2)
    n = (n1 ? 1 : 0) + (n2 ? 1 : 0)
    return SFlat{CAP}(Float64(c), n,
                      ntuple(i -> i == 1 ? (n1 ? Int32(id1) : Int32(id2)) :
                             i == 2 ? (n1 ? Int32(id2) : Int32(0)) : Int32(0), Val(CAP)),
                      ntuple(i -> i == 1 ? (n1 ? Float64(v1) : Float64(v2)) :
                             i == 2 ? (n1 ? Float64(v2) : 0.0) : 0.0, Val(CAP)))
end

function Base.show(io::IO, x::SFlat)
    print(io, "SFlat(", x.c)
    @inbounds for j in 1:x.n
        print(io, " + ", x.vs[j], "ε_", x.ids[j])
    end
    print(io, "; CAP=", length(x.ids), ")")
end
