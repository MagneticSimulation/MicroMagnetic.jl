export build_matrix, build_demag_matrix, is_long_range, LLGJacOperator

import SparseArrays
import LinearAlgebra

# ---- Interaction classification for Route Y -------------------------------
#   Local  → symbolic Epsilon kernel.
#   Demag  → Float64 impulse-columns via effective_field (backend-agnostic).
is_long_range(interaction) = false
is_long_range(interaction::Demag)       = true
is_long_range(interaction::DirectDemag) = true
is_long_range(interaction::DemagFE)     = true


# =============================================================================
# Matrix-free LLG Jacobian operator
# =============================================================================

"""
    LLGJacOperator <: AbstractMatrix{Float64}

Matrix-free linearised LLG Jacobian: dm̄/dt = B · m̄.

Created with `build_matrix(sim; matrixfree=true)`.  Supports:

- `size(op)`, `eltype(op)`, `length(op)`
- `LinearAlgebra.mul!(out, op, x)`  → apply B to any 2N-length x in O(N) space
- `Matrix(op)`, `SparseArrays.sparse(op)`  → materialise (via 2N mul! calls)
- `op[i, j]`  → getindex (via column-j materialisation, O(N) per call, slow)

Implementation (same as Phase-2 Route-Y matrix-free plan):

1. Unwrap tangent x → 3D Cartesian `dm = Σ (x_{2i-1} t_{i,1} + x_{2i} t_{i,2})`
2. Local δH_loc  : single-ε-tag symbolic pass over `m0 + ε·dm` → extract coef
3. Long-range δH_demag : one Float64 `effective_field(demag, sim, dm)` per demag
4. LLG linearization: `δṁ = -γ/(1+α²) · (m0×(δH+α m0×δH) + dm×(H0+α m0×H0))`
5. Project δṁ back onto 2D tangent frame → out vector

So `mul!` uses **O(N)** memory and **1 FFT demag call + 1 local symbolic kernel**
(vs dense Jacobian's 2N calls).
"""
mutable struct LLGJacOperator{T<:AbstractFloat} <: AbstractMatrix{T}
    sim
    N::Int
    gamma::T
    alpha::T
    Rs::Vector{Matrix{T}}          # forward rotation (tangent → Cartesian)
    R_invs::Vector{Matrix{T}}      # inverse rotation (Cartesian → tangent, rows are R^T)
    m0_F64::Vector{T}
    pre_damp::T
    use_damping::Bool
    # Local total H(m0) baseline field
    H0_local::Vector{T}
    # Damping cross-baseline: m0 × H0_local (nothing iff undamped)
    m0xH0_local::Union{Nothing,Vector{T}}
    # Per-demag H0(m0) + m0×H0 if damped (length = number of long-range interactions)
    demag_data::Vector{@NamedTuple{demag, H0::Vector{T}, m0xH0::Union{Nothing,Vector{T}}}}
    # Caches for mul! (avoid per-call 3N allocation on hot Arpack path)
    _dm::Vector{T}
    _δH::Vector{T}
    _w::Vector{T}
    _xcol::Vector{T}                # temporary for Matrix() / getindex() baselines
    _ycol::Vector{T}
end

Base.eltype(::LLGJacOperator{T}) where {T} = T
Base.length(op::LLGJacOperator) = (2 * op.N) ^ 2
Base.size(op::LLGJacOperator) = (2 * op.N, 2 * op.N)
function Base.size(op::LLGJacOperator, d::Integer)
    d in (1, 2) ? 2 * op.N : error("dimension d=$d out of range for size(A, d)")
end

function Base.show(io::IO, ::MIME"text/plain", op::LLGJacOperator)
    println(io, "LLGJacOperator{$(eltype(op))} ($(size(op,1))×$(size(op,2)))")
    println(io, "  gamma       = ", op.gamma)
    println(io, "  alpha       = ", op.alpha)
    println(io, "  N (sites)   = ", op.N)
    println(io, "  damping     = ", op.use_damping ? "Gilbert (α=$(op.alpha))" : "undamped")
    println(io, "  long-range  = ", length(op.demag_data), " interaction(s)")
    print(io,   "  storage     = matrix-free O(N) buffers")
end
Base.show(io::IO, op::LLGJacOperator) = show(io, MIME("text/plain"), op)

# Fallback getindex: materialise the j-th column via mul! on e_j (correct but slow).
function Base.getindex(op::LLGJacOperator{T}, i::Integer, j::Integer) where T
    M, N2 = size(op)
    (1 ≤ i ≤ M && 1 ≤ j ≤ N2) || throw(BoundsError(op, (i, j)))
    e = zeros(T, N2); e[j] = one(T)
    col = similar(e)
    LinearAlgebra.mul!(col, op, e)
    return col[i]
end

# Materialise full 2N×2N dense Jacobian.
function Base.Matrix(op::LLGJacOperator{T}) where T
    M, N2 = size(op)
    B = zeros(T, M, N2)
    e = zeros(T, N2)
    col = similar(e)
    @inbounds for j in 1:N2
        fill!(e, zero(T))
        e[j] = one(T)
        LinearAlgebra.mul!(col, op, e)
        B[:, j] .= col
    end
    return B
end

# Sparse fallback: materialise dense then convert.
function SparseArrays.sparse(op::LLGJacOperator)
    SparseArrays.sparse(Matrix(op))
end

# ----- Linear algebra: mul! --------------------------------------------------

# Internal helper: apply single-tag symbolic local Jacobian `(∂H_loc/∂m)·dm`
function _δH_loc_single!(δH_loc::Vector{Float64}, op::LLGJacOperator, dm::Vector{Float64})
    N = op.N
    m0 = op.m0_F64
    # m = m0 + ε₁·dm  (single tag id=1, unit coefficient 1.0)
    # Note: spin must be AbstractFloat so the KA kernel does T(0) / convert(T, _).
    spin = zeros(AbstractFloat, 3N)
    @inbounds for k in 1:3N
        m0k = m0[k]
        dmk = dm[k]
        spin[k] = ifelse(iszero(dmk), m0k,
                         ifelse(iszero(m0k), Epsilon(1, dmk),
                                m0k + Epsilon(1, dmk)))
    end

    effective_field_local(op.sim, spin)

    # Read ε₁ coefficient of sim.field → δH_loc.
    SymT = Union{Float64, Epsilon, AddExpr}
    @inbounds for k in 1:3N
        fld = convert(SymT, op.sim.field[k])
        δH_loc[k] = get(collect_terms(fld), 1, 0.0)::Float64
    end
    nothing
end

"""
    LinearAlgebra.mul!(out::Vector{<:Real}, op::LLGJacOperator, x::Vector{<:Real})

Apply matrix-free Jacobian `B` to tangent vector `x`, write `out = B · x`.

O(N) memory (all buffers are pre-allocated inside `op`); calls 1 local symbolic
kernel + 1 effective_field per long-range interaction.
"""
function LinearAlgebra.mul!(out::Vector{Tout},
                            op::LLGJacOperator{T},
                            x::Vector{Tx}) where {T<:AbstractFloat,
                                                  Tout<:AbstractFloat,
                                                  Tx<:AbstractFloat}
    N = op.N
    N2 = 2N
    length(x) == N2 && length(out) == N2 ||
        throw(DimensionMismatch("length(out)=$(length(out)), length(x)=$(length(x)), expected $N2"))

    dm = op._dm
    δH = op._δH
    w  = op._w

    # --- 1. Unwrap x: tangent → 3D Cartesian dm -----------------------------
    @inbounds for i in 1:N
        xi = 3(i-1) + 1
        R = op.Rs[i]
        x1 = x[2i-1]; x2 = x[2i]
        dm[xi]   = R[1,1]*x1 + R[1,2]*x2
        dm[xi+1] = R[2,1]*x1 + R[2,2]*x2
        dm[xi+2] = R[3,1]*x1 + R[3,2]*x2
    end

    # --- 2. Local δH_loc via single-ε symbolic pass --------------------------
    _δH_loc_single!(δH, op, dm)

    @inbounds for i in 1:N
        xi = 3(i-1) + 1
        m0x, m0y, m0z = op.m0_F64[xi], op.m0_F64[xi+1], op.m0_F64[xi+2]
        dhx, dhy, dhz = δH[xi],   δH[xi+1],   δH[xi+2]

        px = cross_x(m0x, m0y, m0z, dhx, dhy, dhz)
        py = cross_y(m0x, m0y, m0z, dhx, dhy, dhz)
        pz = cross_z(m0x, m0y, m0z, dhx, dhy, dhz)

        if op.use_damping
            qx = cross_x(m0x, m0y, m0z, px, py, pz)
            qy = cross_y(m0x, m0y, m0z, px, py, pz)
            qz = cross_z(m0x, m0y, m0z, px, py, pz)
            v1 = px + op.alpha * qx
            v2 = py + op.alpha * qy
            v3 = pz + op.alpha * qz
        else
            v1, v2, v3 = px, py, pz
        end

        # dm × H0_local  (+ damping cross terms α·m0×(dm×H0) + α·dm×(m0×H0)
        #                to mirror the Gilbert-damped linearisation used by dense).
        if ms_idx_any_nonzero(op.H0_local, xi)
            hx, hy, hz = op.H0_local[xi], op.H0_local[xi+1], op.H0_local[xi+2]
            dmx, dmy, dmz = dm[xi], dm[xi+1], dm[xi+2]
            v1 += cross_x(dmx, dmy, dmz, hx, hy, hz)
            v2 += cross_y(dmx, dmy, dmz, hx, hy, hz)
            v3 += cross_z(dmx, dmy, dmz, hx, hy, hz)
            if op.use_damping
                # α · m0 × (dm × H0)
                cx = cross_x(dmx, dmy, dmz, hx, hy, hz)
                cy = cross_y(dmx, dmy, dmz, hx, hy, hz)
                cz = cross_z(dmx, dmy, dmz, hx, hy, hz)
                ax = cross_x(m0x, m0y, m0z, cx, cy, cz)
                ay = cross_y(m0x, m0y, m0z, cx, cy, cz)
                az = cross_z(m0x, m0y, m0z, cx, cy, cz)
                v1 += op.alpha * ax
                v2 += op.alpha * ay
                v3 += op.alpha * az
                # α · dm × (m0 × H0)  — cached in m0xH0_local
                m0xH0 = op.m0xH0_local::Vector{T}
                if ms_idx_any_nonzero(m0xH0, xi)
                    bx, by, bz = m0xH0[xi], m0xH0[xi+1], m0xH0[xi+2]
                    v1 += op.alpha * cross_x(dmx, dmy, dmz, bx, by, bz)
                    v2 += op.alpha * cross_y(dmx, dmy, dmz, bx, by, bz)
                    v3 += op.alpha * cross_z(dmx, dmy, dmz, bx, by, bz)
                end
            end
        end
        w[xi], w[xi+1], w[xi+2] = v1, v2, v3
    end

    # --- 3. Long-range: one Float64 effective_field per demag, accumulate ----
    for dd in op.demag_data
        demag, H0_demag, m0xH0 = dd
        fill!(δH, 0.0)
        effective_field(demag, op.sim, dm, 0.0; output=δH)

        @inbounds for i in 1:N
            xi = 3(i-1) + 1
            m0x, m0y, m0z = op.m0_F64[xi], op.m0_F64[xi+1], op.m0_F64[xi+2]
            dhx, dhy, dhz = δH[xi], δH[xi+1], δH[xi+2]

            px = cross_x(m0x, m0y, m0z, dhx, dhy, dhz)
            py = cross_y(m0x, m0y, m0z, dhx, dhy, dhz)
            pz = cross_z(m0x, m0y, m0z, dhx, dhy, dhz)
            if op.use_damping
                qx = cross_x(m0x, m0y, m0z, px, py, pz)
                qy = cross_y(m0x, m0y, m0z, px, py, pz)
                qz = cross_z(m0x, m0y, m0z, px, py, pz)
                v1 = px + op.alpha * qx
                v2 = py + op.alpha * qy
                v3 = pz + op.alpha * qz
            else
                v1, v2, v3 = px, py, pz
            end

            # dm × H0  (+ same Gilbert damping cross terms as local baseline).
            if ms_idx_any_nonzero(H0_demag, xi)
                hx, hy, hz = H0_demag[xi], H0_demag[xi+1], H0_demag[xi+2]
                dmx, dmy, dmz = dm[xi], dm[xi+1], dm[xi+2]
                v1 += cross_x(dmx, dmy, dmz, hx, hy, hz)
                v2 += cross_y(dmx, dmy, dmz, hx, hy, hz)
                v3 += cross_z(dmx, dmy, dmz, hx, hy, hz)
                if op.use_damping
                    # α · m0 × (dm × H0)
                    cx = cross_x(dmx, dmy, dmz, hx, hy, hz)
                    cy = cross_y(dmx, dmy, dmz, hx, hy, hz)
                    cz = cross_z(dmx, dmy, dmz, hx, hy, hz)
                    ax = cross_x(m0x, m0y, m0z, cx, cy, cz)
                    ay = cross_y(m0x, m0y, m0z, cx, cy, cz)
                    az = cross_z(m0x, m0y, m0z, cx, cy, cz)
                    v1 += op.alpha * ax
                    v2 += op.alpha * ay
                    v3 += op.alpha * az
                end
            end
            # α · dm × (m0 × H0)  — cached per-demag (raw cross product; the α factor
            # is multiplied explicitly below together with cross(dm, m0×H0))
            if op.use_damping && m0xH0 !== nothing && ms_idx_any_nonzero(m0xH0, xi)
                ax, ay, az = m0xH0[xi], m0xH0[xi+1], m0xH0[xi+2]
                dmx, dmy, dmz = dm[xi], dm[xi+1], dm[xi+2]
                v1 += op.alpha * cross_x(dmx, dmy, dmz, ax, ay, az)
                v2 += op.alpha * cross_y(dmx, dmy, dmz, ax, ay, az)
                v3 += op.alpha * cross_z(dmx, dmy, dmz, ax, ay, az)
            end

            w[xi]   += v1
            w[xi+1] += v2
            w[xi+2] += v3
        end
    end

    # --- 4. Scale by -γ/(1+α²) & project onto tangent frame ------------------
    gamma_eff = -op.gamma * op.pre_damp
    @inbounds for i in 1:N
        xi = 3(i-1) + 1
        R_inv = op.R_invs[i]
        w1, w2, w3 = w[xi], w[xi+1], w[xi+2]
        out[2i-1] = gamma_eff * (R_inv[1,1]*w1 + R_inv[1,2]*w2 + R_inv[1,3]*w3)
        out[2i]   = gamma_eff * (R_inv[2,1]*w1 + R_inv[2,2]*w2 + R_inv[2,3]*w3)
    end
    return out
end

# 5-argument mul!(y, A, x, α, β) → y = α·A·x + β·y is inherited from the
# AbstractMatrix fallback in LinearAlgebra (calls the 3-arg mul! above and
# fuses the linear combination).  No override needed; saves dispatch ambiguity
# surface with Arpack / KrylovKit wrapper layers.

# --- Arpack / Krylov traits --------------------------------------------------

# LLG Jacobian is NOT symmetric in general. Explicit false to avoid solvers
# taking a symmetric shortcut that would break for damped / non-square systems.
LinearAlgebra.issymmetric(::LLGJacOperator) = false
LinearAlgebra.ishermitian(op::LLGJacOperator) = issymmetric(op)
LinearAlgebra.issuccess(::LLGJacOperator) = true

"""
    copy(op::LLGJacOperator) -> LLGJacOperator

Shallow copy of operator metadata with freshly allocated caches so that
multiple Arpack solves (or threaded multiplies) do not race on `_dm/_δH/…`.
Heavy read-only objects (sim, demag plans, rotation matrices, baselines)
are shared between copies.
"""
function Base.copy(op::LLGJacOperator{T}) where T
    n3 = 3 * op.N
    n2 = 2 * op.N
    LLGJacOperator{T}(
        op.sim, op.N, op.gamma, op.alpha,
        op.Rs, op.R_invs,          # rotation matrices: immutable after build → share
        op.m0_F64,                 # baseline spin: read-only → share
        op.pre_damp, op.use_damping,
        op.H0_local,               # read-only baseline → share
        op.m0xH0_local,            # damping cross-baseline: read-only → share
        op.demag_data,             # (demag ref, H0 vec, m0xH0 vec): read-only → share
        zeros(T, n3), zeros(T, n3), zeros(T, n3),  # _dm, _δH, _w: FRESH
        zeros(T, n2), zeros(T, n2))                 # _xcol, _ycol: FRESH
end


# Hot helper: "any of 3 consecutive entries nonzero" — used to skip the zero
# baseline-field cross product branches on sites that carry no contribution.
@inline function ms_idx_any_nonzero(vec::Vector{T}, i3::Int) where T<:AbstractFloat
    @inbounds return !(iszero(vec[i3]) && iszero(vec[i3+1]) && iszero(vec[i3+2]))
end

# ----- Construct operator from sim ------------------------------------------

function _make_operator(sim; gamma, alpha)
    N = sim.n_total
    T = Float64
    use_damping = abs(alpha) > eps(T)
    @info("Building LLGJacOperator (matrix-free, O(N) storage)  " *
          "N=$N  gamma=$gamma  alpha=$alpha  damping=$use_damping")
    if eltype(sim.spin) !== AbstractFloat
        @warn("""Operator mode requires MicroMagnetic.set_precision(AbstractFloat)
                 BEFORE constructing the sim; eltype(sim.spin) = $(eltype(sim.spin)).
                 Without it the single-ε local symbolic kernel may fail.""")
    end

    m0_F64::Vector{T} = convert(Vector{T}, Array(sim.spin))

    Rs     = Vector{Matrix{T}}(undef, N)
    R_invs = similar(Rs)
    @inbounds for i in 1:N
        x = 3(i-1) + 1
        R = rotation_matrix(m0_F64[x], m0_F64[x+1], m0_F64[x+2])
        Rs[i]     = R
        R_invs[i] = [R[1,1] R[2,1] R[3,1];
                     R[1,2] R[2,2] R[3,2];
                     R[1,3] R[2,3] R[3,3]]
    end

    # Local + Zeeman baseline at m0 (Float64 spin → baseline field)
    H0_local = zeros(T, 3N)
    _compute_baseline_local!(H0_local, sim, m0_F64)

    # Damping baseline for local interactions: compute m0 × H0_local now so
    # mul! can reuse it for the α·dm×(m0×H0) Gilbert cross term.
    m0xH0_local = if use_damping
        out = zeros(T, 3N)
        @inbounds for i in 1:N
            xi = 3(i-1) + 1
            a, b, c = m0_F64[xi], m0_F64[xi+1], m0_F64[xi+2]
            d, e, f = H0_local[xi], H0_local[xi+1], H0_local[xi+2]
            out[xi]   = cross_x(a, b, c, d, e, f)
            out[xi+1] = cross_y(a, b, c, d, e, f)
            out[xi+2] = cross_z(a, b, c, d, e, f)
        end
        out
    else
        nothing
    end

    # Long-range baseline H0 + m0×H0 (damped) per demag
    lr_list = [inter for inter in sim.interactions if is_long_range(inter)]
    demag_data = map(lr_list) do demag
        H0 = zeros(T, 3N)
        effective_field(demag, sim, m0_F64, 0.0; output=H0)
        m0xH0 = if use_damping
            out = zeros(T, 3N)
            @inbounds for i in 1:N
                xi = 3(i-1)+1
                a, b, c = m0_F64[xi], m0_F64[xi+1], m0_F64[xi+2]
                d, e, f = H0[xi],    H0[xi+1],    H0[xi+2]
                out[xi]   = cross_x(a, b, c, d, e, f)
                out[xi+1] = cross_y(a, b, c, d, e, f)
                out[xi+2] = cross_z(a, b, c, d, e, f)
            end
            out
        else
            nothing
        end
        (; demag, H0, m0xH0)
    end

    pre_damp = 1.0 / (1.0 + alpha*alpha)
    n3 = 3N
    n2 = 2N
    _dm   = zeros(T, n3)
    _δH   = zeros(T, n3)
    _w    = zeros(T, n3)
    _xcol = zeros(T, n2)
    _ycol = zeros(T, n2)
    return LLGJacOperator{T}(sim, N, gamma, alpha, Rs, R_invs, m0_F64,
                             pre_damp, use_damping,
                             H0_local, m0xH0_local, demag_data,
                             _dm, _δH, _w, _xcol, _ycol)
end

# Float64 baseline H0 of local + Zeeman interactions evaluated at m0_F64.
function _compute_baseline_local!(out::Vector{Float64}, sim, m0_F64::Vector{Float64})
    N = length(m0_F64) ÷ 3
    local_spin = Vector{AbstractFloat}(undef, 3N)
    @inbounds for k in 1:(3N); local_spin[k] = m0_F64[k]; end
    effective_field_local(sim, local_spin)
    # Use cpart() rather than convert(Vector{Float64}, …): the latter would
    # call Float64(::FlatTerm/Epsilon) per element and trigger the loud
    # convert stub in src/eigen/util.jl.  In the baseline path field[k] is
    # always Float64-typed at runtime, but cpart is the safe, future-proof
    # accessor that documents intent ("extract the ε⁰ constant part").
    copyto!(out, cpart(sim.field))
    nothing
end


# =============================================================================
# build_matrix public API — now with matrixfree kwarg
# =============================================================================

@doc raw"""
    build_matrix(sim; gamma=2.21e5, sparse=false, alpha=0.01, matrixfree=false)

Build the linearised LLG Jacobian `B` that maps tangent perturbations
`dm̄` to `dm̄/dt = B · dm̄`.

# Precision prerequisite

Local interactions use symbolic `Epsilon` / `FlatTerm` duals (see
`eigen/util.jl`).  The pre-allocated field buffers inside `sim` and its
interactions must therefore be `AbstractFloat` typed, i.e. the sim must be
constructed after calling `MicroMagnetic.set_precision(AbstractFloat)`.

# Output modes

| kwarg              | type returned                       | memory   | cost per call                  |
| :----------------- | :---------------------------------- | :------- | :----------------------------- |
| (default)          | `Matrix{Float64}`                   | O(4N²)   | 1 sym-pass + 2N demag calls    |
| `sparse=true`      | `SparseMatrixCSC{Float64,Int}`      | ~O(NZ)   | as dense then convert          |
| `matrixfree=true`  | `LLGJacOperator` (AbstractMatrix)   | O(N)     | `mul!`: 1 sym-pass + 1 demag   |

# Large-scale eigenmode calculation with Arpack

For meshes where the dense 2N×2N Jacobian does not fit in memory, pass
`matrixfree=true` and use `Arpack.eigs` on the returned operator.  The
operator fully implements the `AbstractMatrix` interface Arpack needs
(`size`, `eltype`, 3- and 5-argument `mul!`, `issymmetric`, `copy`).

```julia
using Arpack
using MicroMagnetic

MicroMagnetic.set_precision(AbstractFloat)

mesh = FDMesh(nx=8, ny=8, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
sim  = Sim(mesh); set_Ms(sim, 8e5); init_m0(sim, (1,0,0))
add_exch(sim, 1.3e-11); add_demag(sim); add_zeeman(sim, (0,0,1000))

op = build_matrix(sim; matrixfree=true, alpha=0.01)

# Solve B·v = λ·v for the 20 eigenvalues of largest real part
# (→ fastest-growing / least-damped eigenmodes):
vals, vecs, nconv, niter, nmult, resid =
    Arpack.eigs(op; nev=20, which=:LR, tol=1e-8, maxiter=1000,
                ritzvec=true, explicittransform=:none)

# Convert angular eigenvalue σ (damping, rad/s) to frequency f (GHz)
# and relaxation rate Γ: f = |imag(σ)|/(2π·1e9), Γ = -real(σ)/(2π·1e9).
freq_GHz    = abs.(imag.(vals)) ./ (2π * 1e9)
damping_GHz = -real.(vals) ./ (2π * 1e9)
```

For undamped eigenproblems (`alpha=0.0`) `B` is a real linearised Landau-Lifshitz
operator and you may want `which=:LM` instead of `:LR` to obtain the highest
swing frequencies first.  For non-Hermitian / damped cases `:LR` (largest real
part) is the standard choice to pick physically interesting modes.
"""
function build_matrix(sim; gamma=2.21e5, sparse=false, alpha=0.01, matrixfree=false)
    if eltype(sim.spin) !== AbstractFloat
        @warn("""build_matrix requires MicroMagnetic.set_precision(AbstractFloat)
                 BEFORE constructing the sim; current eltype(sim.spin) = $(eltype(sim.spin)).
                 Otherwise symbolic writes fail with MethodError: Float64(::AddExpr).""")
    end
    if matrixfree
        sparse && @info("  sparse=true ignored with matrixfree=true; use " *
                        "SparseArrays.sparse(op) or Matrix(op) on demand.")
        return _make_operator(sim; gamma=gamma, alpha=alpha)
    end
    use_dm_dense = abs(alpha) > eps(Float64)
    @info("Building dense LLG Jacobian (Route Y: Epsilon local + Float64 demag, 2N×2N)  " *
          "N=$(sim.n_total)  gamma=$gamma  alpha=$alpha  damping=$use_dm_dense")
    return _build_matrix_impl(sim; gamma=gamma, sparse=sparse, alpha=alpha)
end


function _build_matrix_impl(sim; gamma, sparse, alpha)
    N = sim.n_total
    T = Float64
    pre_damp    = 1.0 / (1.0 + alpha*alpha)
    use_damping = abs(alpha) > eps(T)
    @info("  _build_matrix_impl:  received gamma=$gamma  alpha=$alpha  " *
          "pre_damp=$pre_damp  use_damping=$use_damping")
    local SymT = Union{Float64, Epsilon, AddExpr}

    m0_F64::Vector{Float64} = cpart(sim.spin)
    Ms_raw = isa(sim, MicroSim) ? Array(sim.mu0_Ms)::Vector{<:AbstractFloat} :
                                  Array(sim.mu_s)::Vector{<:AbstractFloat}
    Ms_F64::Vector{Float64} = cpart(Ms_raw)

    Rs     = Vector{Matrix{Float64}}(undef, N)
    R_invs = similar(Rs)
    @inbounds for i in 1:N
        x = 3(i-1) + 1
        R = rotation_matrix(m0_F64[x], m0_F64[x+1], m0_F64[x+2])
        Rs[i]     = R
        R_invs[i] = [R[1,1] R[2,1] R[3,1];
                     R[1,2] R[2,2] R[3,2];
                     R[1,3] R[2,3] R[3,3]]
    end

    # Step 1: spin = m0 + Σ_j (ε_{2j-1} t_{j,1} + ε_{2j} t_{j,2})
    spin  = zeros(AbstractFloat, 3N)
    dm_dt = Vector{SymT}(undef, 3N)
    fill!(dm_dt, 0.0)

    for i in 1:N
        if @inbounds Ms_F64[i] == 0
            @inbounds(spin[3i-2] = 0.0); @inbounds(spin[3i-1] = 0.0); @inbounds(spin[3i] = 0.0)
            continue
        end
        x = 3(i-1)
        R = Rs[i]
        v1 = Epsilon(2i-1, 1.0)
        v2 = Epsilon(2i,   1.0)
        @inbounds begin
            spin[x+1] = R[1,1]*v1 + R[1,2]*v2 + R[1,3]
            spin[x+2] = R[2,1]*v1 + R[2,2]*v2 + R[2,3]
            spin[x+3] = R[3,1]*v1 + R[3,2]*v2 + R[3,3]
        end
    end
    for k in eachindex(spin); spin[k] = simplify(spin[k]); end

    # Step 2a: local → symbolic H_eff
    @info("  Local: 2N-ε symbolic pass (Exchange/Anisotropy/DMI/Zeeman …) ...")
    effective_field_local(sim, spin)
    H_eff = Vector{SymT}(undef, 3N)
    @inbounds for k in eachindex(H_eff)
        H_eff[k] = simplify(convert(SymT, sim.field[k]))
    end

    # Step 3: dmdt = -γ (m × H + α m × (m × H))   (pre_damp applied at end)
    @info("  Assembling m × H_eff (symbolic) ...")
    for i in 1:N
        x = 3(i-1)
        @inbounds sx, sy, sz = spin[x+1], spin[x+2], spin[x+3]
        @inbounds hx, hy, hz = H_eff[x+1], H_eff[x+2], H_eff[x+3]
        px = cross_x(sx, sy, sz, hx, hy, hz)
        py = cross_y(sx, sy, sz, hx, hy, hz)
        pz = cross_z(sx, sy, sz, hx, hy, hz)
        if use_damping
            qx = cross_x(sx, sy, sz, px, py, pz)
            qy = cross_y(sx, sy, sz, px, py, pz)
            qz = cross_z(sx, sy, sz, px, py, pz)
            dm_dt[x+1] = -gamma * (px + alpha * qx)
            dm_dt[x+2] = -gamma * (py + alpha * qy)
            dm_dt[x+3] = -gamma * (pz + alpha * qz)
        else
            dm_dt[x+1] = -gamma * px
            dm_dt[x+2] = -gamma * py
            dm_dt[x+3] = -gamma * pz
        end
    end
    for k in eachindex(dm_dt); dm_dt[k] = simplify(dm_dt[k]); end

    # Step 4: project onto tangent frame
    for i in 1:N
        x = 3(i-1)
        R_inv = R_invs[i]
        @inbounds v1, v2, v3 = dm_dt[x+1], dm_dt[x+2], dm_dt[x+3]
        @inbounds begin
            spin[x+1] = R_inv[1,1]*v1 + R_inv[1,2]*v2 + R_inv[1,3]*v3
            spin[x+2] = R_inv[2,1]*v1 + R_inv[2,2]*v2 + R_inv[2,3]*v3
            spin[x+3] = R_inv[3,1]*v1 + R_inv[3,2]*v2 + R_inv[3,3]*v3
        end
    end
    for k in eachindex(spin); spin[k] = simplify(spin[k]); end

    # Step 5: ε-coefficients → dense B_loc
    @info("  Extracting local ε-coefficients → 2N×2N block ...")
    B = zeros(Float64, 2N, 2N)
    for i in 1:N
        x = 3(i-1)
        terms1 = collect_terms(spin[x+1])
        @inbounds for j in 1:(2N)
            v = get(terms1, j, 0.0)::Float64
            abs(v) > eps(Float64) && (B[2i-1, j] = v)
        end
        terms2 = collect_terms(spin[x+2])
        @inbounds for j in 1:(2N)
            v = get(terms2, j, 0.0)::Float64
            abs(v) > eps(Float64) && (B[2i, j] = v)
        end
    end
    use_damping && (B .*= pre_damp)

    # Step 2b: demag — Float64 impulse columns
    lr_list = [inter for inter in sim.interactions if is_long_range(inter)]
    if !isempty(lr_list)
        @info("  Demag: Float64 impulse-columns (2N × effective_field, backend-agnostic) ...")
        for demag in lr_list
            add_demag_block!(B, demag, sim, m0_F64, Ms_F64, Rs, R_invs;
                             gamma=gamma, alpha=alpha, pre_damp=pre_damp,
                             use_damping=use_damping)
        end
    end

    sparse && @info("  Converting dense result to SparseMatrixCSC ...")
    return sparse ? SparseArrays.sparse(B) : B
end


"""
Accumulate effective fields of LOCAL interactions (skip `is_long_range`) into
`sim.field`.  Requires sim created with `set_precision(AbstractFloat)`.
"""
function effective_field_local(sim, spin)
    fill!(sim.field, 0.0)
    for interaction in sim.interactions
        is_long_range(interaction) && continue
        effective_field(interaction, sim, spin, 0.0)
        vector_add(sim.field, interaction.field)
    end
    nothing
end


"""
In-place add the demag contribution of a single `Demag` backend to `matrix`.
"""
function add_demag_block!(matrix::Matrix{Float64}, demag, sim,
                          m0_F64::Vector{Float64}, Ms_F64::Vector{Float64},
                          Rs::Vector{Matrix{Float64}}, R_invs::Vector{Matrix{Float64}};
                          gamma=2.21e5, alpha=0.0, pre_damp=1.0, use_damping=false)
    N = sim.n_total

    H0 = zeros(Float64, 3N)
    effective_field(demag, sim, m0_F64, 0.0; output=H0)

    m0xH0 = if use_damping
        out = zeros(Float64, 3N)
        @inbounds for i in 1:N
            xi = 3(i-1) + 1
            a, b, c = m0_F64[xi], m0_F64[xi+1], m0_F64[xi+2]
            d, e, f = H0[xi],   H0[xi+1],   H0[xi+2]
            out[xi]   = cross_x(a, b, c, d, e, f)
            out[xi+1] = cross_y(a, b, c, d, e, f)
            out[xi+2] = cross_z(a, b, c, d, e, f)
        end
        out
    else
        nothing
    end

    dm = zeros(Float64, 3N)
    dH = zeros(Float64, 3N)
    gamma_eff = -gamma * pre_damp
    alpha_eff = alpha

    @inbounds for j in 1:N
        Ms_F64[j] == 0 && continue
        xj = 3(j-1) + 1
        Rj = Rs[j]
        for alpha_axis in 1:2
            c = 2(j-1) + alpha_axis
            t1u = Rj[1, alpha_axis]; t2u = Rj[2, alpha_axis]; t3u = Rj[3, alpha_axis]
            fill!(dm, 0.0)
            dm[xj] = t1u; dm[xj+1] = t2u; dm[xj+2] = t3u

            fill!(dH, 0.0)
            effective_field(demag, sim, dm, 0.0; output=dH)

            for i in 1:N
                xi = 3(i-1) + 1
                m0x, m0y, m0z = m0_F64[xi], m0_F64[xi+1], m0_F64[xi+2]
                dhx, dhy, dhz = dH[xi],   dH[xi+1],   dH[xi+2]

                px = cross_x(m0x, m0y, m0z, dhx, dhy, dhz)
                py = cross_y(m0x, m0y, m0z, dhx, dhy, dhz)
                pz = cross_z(m0x, m0y, m0z, dhx, dhy, dhz)
                if use_damping
                    qx = cross_x(m0x, m0y, m0z, px, py, pz)
                    qy = cross_y(m0x, m0y, m0z, px, py, pz)
                    qz = cross_z(m0x, m0y, m0z, px, py, pz)
                    v1 = px + alpha_eff * qx
                    v2 = py + alpha_eff * qy
                    v3 = pz + alpha_eff * qz
                else
                    v1, v2, v3 = px, py, pz
                end
                if i == j
                    h0x, h0y, h0z = H0[xi], H0[xi+1], H0[xi+2]
                    # dm × H0 — shared between conservative and one damping cross
                    cx = cross_x(t1u, t2u, t3u, h0x, h0y, h0z)
                    cy = cross_y(t1u, t2u, t3u, h0x, h0y, h0z)
                    cz = cross_z(t1u, t2u, t3u, h0x, h0y, h0z)
                    v1 += cx
                    v2 += cy
                    v3 += cz
                    if use_damping
                        # α · m0 × (dm × H0)  — required whenever H0∥m0 (demag self-field
                        # case) because m0×H0=0 then but this term still contributes a
                        # diagonal damping component (mirrors matrix-free mul! lines 263-273)
                        # (m0x/m0y/m0z re-used from outer per-site declaration at L707)
                        ax = cross_x(m0x, m0y, m0z, cx, cy, cz)
                        ay = cross_y(m0x, m0y, m0z, cx, cy, cz)
                        az = cross_z(m0x, m0y, m0z, cx, cy, cz)
                        v1 += alpha_eff * ax
                        v2 += alpha_eff * ay
                        v3 += alpha_eff * az
                        # α · dm × (m0 × H0)  — zero when m0∥H0, non-zero otherwise
                        if m0xH0 !== nothing
                            bx, by, bz = m0xH0[xi], m0xH0[xi+1], m0xH0[xi+2]
                            v1 += alpha_eff * cross_x(t1u, t2u, t3u, bx, by, bz)
                            v2 += alpha_eff * cross_y(t1u, t2u, t3u, bx, by, bz)
                            v3 += alpha_eff * cross_z(t1u, t2u, t3u, bx, by, bz)
                        end
                    end
                end
                R_inv = R_invs[i]
                fx = gamma_eff * (R_inv[1,1]*v1 + R_inv[1,2]*v2 + R_inv[1,3]*v3)
                fy = gamma_eff * (R_inv[2,1]*v1 + R_inv[2,2]*v2 + R_inv[2,3]*v3)
                matrix[2i-1, c] += fx
                matrix[2i,     c] += fy
            end
        end
    end
    nothing
end


"""
    build_demag_matrix(demag, sim; gamma=2.21e5, sparse=false, alpha=0.01)

Demag-only 2N×2N LLG Jacobian.  Pure Float64, zero symbolic overhead;
`sparse=true` returns a `SparseMatrixCSC`.
"""
function build_demag_matrix(demag, sim; gamma=2.21e5, sparse=false, alpha=0.01)
    @info("Building demag-only LLG Jacobian (impulse-columns, dense 2N×2N) ...")
    N = sim.n_total
    m0_F64::Vector{Float64} = cpart(sim.spin)
    Ms_raw = isa(sim, MicroSim) ? Array(sim.mu0_Ms)::Vector{<:AbstractFloat} :
                                  Array(sim.mu_s)::Vector{<:AbstractFloat}
    Ms_F64::Vector{Float64} = cpart(Ms_raw)

    Rs     = [rotation_matrix(m0_F64[3i-2], m0_F64[3i-1], m0_F64[3i]) for i in 1:N]
    R_invs = map(R -> [R[1,1] R[2,1] R[3,1];
                       R[1,2] R[2,2] R[3,2];
                       R[1,3] R[2,3] R[3,3]], Rs)

    pre_damp    = 1.0 / (1.0 + alpha*alpha)
    use_damping = abs(alpha) > eps(Float64)

    B = zeros(Float64, 2N, 2N)
    add_demag_block!(B, demag, sim, m0_F64, Ms_F64, Rs, R_invs;
                     gamma=gamma, alpha=alpha, pre_damp=pre_damp,
                     use_damping=use_damping)

    sparse && @info("  Converting dense result to SparseMatrixCSC ...")
    return sparse ? SparseArrays.sparse(B) : B
end
