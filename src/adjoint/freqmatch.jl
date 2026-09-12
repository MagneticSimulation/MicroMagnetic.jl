# FrequencyMatch scenario — eigenvalue adjoint for FMR frequency matching.
#
# Forward  (run_forward!): dense B = build_matrix(sim; gamma, alpha) at the
#     frozen equilibrium m0, eigenpair (λ, x) selected by |Im λ| closest to the
#     target frequency (the exact ±iω tie of a real B resolves to the +Im
#     branch), left eigenvector y = eigenvector of Bᵀ for the same λ; loss
#     J = (ω_sel − ω*)² with ω = Im(λ)/(2π·1e9) in GHz (test_cubic convention).
# Backward (gradient!):    dλ/dθ = yᵀ(∂B/∂θ)x / (yᵀx)  and
#     dJ/dθ = 2(ω_sel − ω*)·Im(dλ/dθ)/(2π·1e9).  ∂B/∂θ is assembled exactly —
#     B is strictly affine in the scalar Ku/Aex/D design variables of P2a, so
#     two unit-parameter build_matrix calls give ∂B/∂θ = B(θ+1) − B(θ) with no
#     finite-difference truncation error.
#
# Dense path for small N; the non-uniform P5 extension (matrixfree + Arpack +
# LLGJacOperator adjoint) reuses the same dλ/dθ formula.

export FrequencyMatch, run_forward!, gradient!, set_design!, value_and_gradient

mutable struct FrequencyMatch
    sim::AbstractSim
    design::NamedTuple           # θ source of truth: (Ku=…, Aex=…, D=…) scalars
    target_freq::Float64         # GHz
    gamma::Float64
    alpha::Float64
    slots::NamedTuple            # design key → the interaction its writes go to
    # forward cache (all keyed to the current design point; nothing = stale)
    B
    lambda
    xsel
    ysel
    omega::Float64               # selected mode frequency, GHz (signed Im λ)
    J
    FrequencyMatch() = new()     # internal: filled by the constructor below
end

@doc raw"""
    FrequencyMatch(sim; design, target, gamma=2.21e5, alpha=0.0)

Frequency-matching inverse problem: find the scalar design variables θ
that drive a selected eigenmode of the linearised LLG operator `B` (built at
the frozen equilibrium `sim.spin`) to a target frequency.

```julia
fm = FrequencyMatch(sim; design = (Ku = K0,), target = (freq = 8.0,))
J = run_forward!(fm)      # J = (ω_sel − ω*)²,  ω in GHz
G = gradient!(fm)         # NamedTuple shaped like `design`, e.g. G.Ku
set_design!(fm; Ku = K1)  # write a new design point (invalidates the cache)
```

`target` is the target frequency in GHz — a bare Real (`target = 8.0`) or a
NamedTuple with a `:freq` entry.

# Design variables

`design` takes **scalars** from `Ku` (uniaxial anisotropy, `add_anis`),
`Aex` (isotropic exchange, `add_exch`) and `D` (isotropic DMI, `add_dmi`,
Dx=Dy=Dz moving together); each requires exactly one matching interaction in
`sim`.  `B` is strictly affine in these, which is what makes the unit-parameter
`∂B/∂θ = B(θ+1) − B(θ)` assembly exact.  **`Ms` is not supported**: ω depends on
Ms nonlinearly through the field denominators, which needs the quotient rule
(planned for a later version).  The declared values are written into `sim` at
construction and on every [`set_design!`](@ref) / [`run_forward!`](@ref) — the
problem owns its parameters from then on.

# Selection convention

The eigenmode is the eigenvalue with |Im λ| closest to `target.freq·2π·1e9`
(rad/s); the exact ±iω tie of a real `B` resolves to the +Im branch, so
`ω = Im(λ)/(2π·1e9)` (GHz, `test_cubic` convention) is consistently the positive
branch.

# Precision prerequisite

Dense `build_matrix` runs the symbolic local pass, so `sim` must be constructed
**after** `MicroMagnetic.set_precision(AbstractFloat)` (inherited from the
eigen layer; lifted once shadow-sim automation lands).

# Gradient

`dλ/dθ = yᵀ(∂B/∂θ)x / (yᵀx)` with `y` the left eigenvector (eigenvector of Bᵀ
for the same λ) and `dJ/dθ = 2(ω_sel − ω*)·Im(dλ/dθ)/(2π·1e9)` — one backward
pass per gradient, independent of the number of design variables.  Dense `B`
means small N; memory is O(4N²).
"""
function FrequencyMatch(sim::AbstractSim; design, target, gamma=2.21e5, alpha=0.0)
    design isa NamedTuple && !isempty(design) ||
        throw(ArgumentError("design must be a non-empty NamedTuple of scalar design " *
                            "variables, e.g. design=(Ku = 5e4,)"))
    if :Ms in keys(design)
        throw(ArgumentError("Ms is not a supported FrequencyMatch design variable in P2a: " *
                            "ω depends on Ms nonlinearly (demag/anisotropy denominators); " *
                            "the quotient rule lands in a later version."))
    end
    for (k, v) in pairs(design)
        k in (:Ku, :Aex, :D) ||
            throw(ArgumentError("unsupported design variable :$k — P2a supports scalar :Ku, :Aex, :D"))
        v isa Real ||
            throw(ArgumentError("design variable :$k must be a scalar Real in P2a (got $(typeof(v)))"))
    end
    (target isa Real || (target isa NamedTuple && :freq in keys(target))) ||
        throw(ArgumentError("target must be the target frequency in GHz — a bare Real " *
                            "(target = 8.0) or a NamedTuple with a :freq entry " *
                            "(target = (freq = 8.0,))"))

    slots = NamedTuple{keys(design)}(Tuple(_resolve_slot(sim, k) for k in keys(design)))
    fm = FrequencyMatch()
    fm.sim = sim
    fm.design = NamedTuple{keys(design)}(Tuple(Float64(v) for v in values(design)))
    fm.target_freq = Float64(target isa Real ? target : target.freq)
    fm.gamma = Float64(gamma)
    fm.alpha = Float64(alpha)
    fm.slots = slots
    fm.B = nothing; fm.lambda = nothing; fm.xsel = nothing; fm.ysel = nothing
    fm.omega = NaN; fm.J = nothing

    for (k, θ) in pairs(fm.design)          # declared values become the sim's values
        _write_design!(slots[k], θ)
    end
    return fm
end

# ----------------------------------------------------------------------------------
# Design slots: resolve the interaction each design variable writes into
# ----------------------------------------------------------------------------------

function _only_interaction(sim, what::Symbol, pred, hint::String)
    cands = filter(pred, sim.interactions)
    isempty(cands) && throw(ArgumentError(
        "FrequencyMatch: no interaction matches design variable :$what ($hint)."))
    length(cands) > 1 && throw(ArgumentError(
        "FrequencyMatch: design variable :$what is ambiguous — $(length(cands)) matching " *
        "interactions; P2a supports exactly one."))
    return cands[1]
end

function _resolve_slot(sim, kind::Symbol)
    if kind === :Ku
        return _only_interaction(sim, kind, i -> i isa Anisotropy,
                                 "add one uniaxial anisotropy via add_anis(sim, Ku)")
    elseif kind === :Aex
        return _only_interaction(sim, kind, i -> i isa Exchange && i.Ax === i.Ay === i.Az,
                                 "add isotropic exchange via add_exch(sim, A) with a single number")
    elseif kind === :D
        return _only_interaction(sim, kind, i -> i isa DMI && i.Dx === i.Dy === i.Dz,
                                 "add isotropic DMI via add_dmi(sim, D) with a single number")
    end
    throw(ArgumentError("unsupported design variable :$kind"))  # unreachable (validated up front)
end

# Fill params are immutable — replace; dense arrays are filled in place.
_design_array(src::Fill, v::Float64) = Fill(v, length(src))
_design_array(src::AbstractArray, v::Float64) = (fill!(src, v); src)

# In-place parameter writes must drop the stencil-layer pair tables: those key
# on mesh.layout_version only (src/micro/stencil.jl), not on parameter values.
_invalidate_design_caches!(::Anisotropy) = nothing
function _invalidate_design_caches!(inter::Exchange)
    inter.pair_Ax = inter.pair_Ay = inter.pair_Az = nothing
    inter.stencil_layout = -1
end
function _invalidate_design_caches!(inter::DMI)
    inter.pair_Dx = inter.pair_Dy = inter.pair_Dz = nothing
    empty!(inter.Dcls)
    inter.stencil_layout = -1
end

function _write_design!(inter::Anisotropy, v::Float64)
    inter.Ku = _design_array(inter.Ku, v)
    _invalidate_design_caches!(inter)
end
function _write_design!(inter::Exchange, v::Float64)
    inter.Ax = _design_array(inter.Ax, v)
    inter.Ay = inter.Az = inter.Ax          # scalar Aex keeps the three bonds aliased
    _invalidate_design_caches!(inter)
end
function _write_design!(inter::DMI, v::Float64)
    inter.Dx = _design_array(inter.Dx, v)
    inter.Dy = inter.Dz = inter.Dx          # scalar D moves Dx=Dy=Dz together
    _invalidate_design_caches!(inter)
end

# ----------------------------------------------------------------------------------
# Forward: B, eigenpair selection, left eigenvector
# ----------------------------------------------------------------------------------

# |Im λ| closest to the target angular frequency; inside a tie window (the ±iω
# pair of a real B ties exactly) prefer the +Im branch so ω is the positive one.
function _select_mode(evals, target_rad::Float64)
    d = abs.(abs.(imag.(evals)) .- target_rad)
    window = 1e-8 * max(target_rad, 1.0)
    near = findall(x -> x <= minimum(d) + window, d)
    imax = near[1]
    for j in near
        imag(evals[j]) > imag(evals[imax]) && (imax = j)
    end
    return imax
end

"""
    run_forward!(fm::FrequencyMatch) -> J

Sync the declared design variables into `sim`, build the dense linearised
operator `B = build_matrix(sim; gamma, alpha)`, select the eigenmode by |Im λ|
closest to the target frequency (ties → +Im branch), and cache `(λ, x, y, B)`;
the left eigenvector `y` is the eigenvector of `Bᵀ` for the same λ.  Returns
`J = (ω_sel − ω*)²` with `ω_sel = Im(λ)/(2π·1e9)` in GHz.
"""
function run_forward!(fm::FrequencyMatch)
    for (k, θ) in pairs(fm.design)
        _write_design!(fm.slots[k], θ)
    end
    B = build_matrix(fm.sim; gamma = fm.gamma, alpha = fm.alpha)
    EB = eigen(B)
    target_rad = 2π * 1e9 * fm.target_freq
    j = _select_mode(EB.values, target_rad)
    lambda = EB.values[j]
    x = EB.vectors[:, j]

    BT = Matrix(transpose(B))
    ET = eigen(BT)
    k = argmin(abs.(ET.values .- lambda))
    abs(ET.values[k] - lambda) <= 1e-6 * (1 + abs(lambda)) ||
        throw(ErrorException("FrequencyMatch: left-eigenvector pairing failed " *
            "(|λ_Bᵀ − λ_B| = $(abs(ET.values[k] - lambda))); B looks defective."))
    y = ET.vectors[:, k]
    abs(transpose(y) * x) > 1e-10 || throw(ErrorException(
        "FrequencyMatch: selected eigenvalue λ = $lambda has yᵀx ≈ 0 (defective or " *
        "degenerate mode); the eigenvalue adjoint is ill-defined there."))

    omega = imag(lambda) / (2π * 1e9)
    omega == 0.0 && fm.target_freq != 0.0 && @warn(
        "FrequencyMatch: the selected eigenvalue is real (|Im λ| = 0) — no precessional " *
        "mode matches the target; check alpha/gamma and the target frequency.")

    fm.B = B; fm.lambda = lambda; fm.xsel = x; fm.ysel = y
    fm.omega = omega
    fm.J = (omega - fm.target_freq)^2
    return fm.J
end

# ----------------------------------------------------------------------------------
# Backward: eigenvalue adjoint with exact unit-parameter ∂B/∂θ
# ----------------------------------------------------------------------------------

"""
    gradient!(fm::FrequencyMatch) -> NamedTuple

Eigenvalue adjoint at the cached forward point (runs [`run_forward!`](@ref)
first if the cache is empty).  For each design variable θ:

```
∂B/∂θ = B(θ+1) − B(θ)                     # exact: B is affine in scalar Ku/Aex/D
dλ/dθ = yᵀ(∂B/∂θ)x / (yᵀx)                # y = left eigenvector from Bᵀ
dJ/dθ = 2(ω_sel − ω*)·Im(dλ/dθ)/(2π·1e9)  # ω = Im(λ)/(2π·1e9) GHz
```

Returns a NamedTuple shaped like `design`.  The `sim` parameters are restored
to the design point afterwards.
"""
function gradient!(fm::FrequencyMatch)
    fm.J === nothing && run_forward!(fm)
    B, x, y = fm.B, fm.xsel, fm.ysel
    omega_err = fm.omega - fm.target_freq
    G = Float64[]
    for (k, θ0) in pairs(fm.design)
        _write_design!(fm.slots[k], θ0 + 1.0)               # unit parameter
        Bp = build_matrix(fm.sim; gamma = fm.gamma, alpha = fm.alpha)
        dlambda = transpose(y) * ((Bp - B) * x) / (transpose(y) * x)
        domega = imag(dlambda) / (2π * 1e9)
        push!(G, 2 * omega_err * domega)
        _write_design!(fm.slots[k], θ0)                     # restore
    end
    return NamedTuple{keys(fm.design)}(Tuple(G))
end

"""
    set_design!(fm::FrequencyMatch; kw...) -> design

Write new scalar values for declared design variables (e.g.
`set_design!(fm; Ku = 6e4)`), sync them into `sim`, and invalidate the forward
cache — the next [`run_forward!`](@ref) / [`gradient!`](@ref) rebuilds at the
new point.  Unknown keys are rejected.
"""
function set_design!(fm::FrequencyMatch; kw...)
    for (k, v) in pairs(kw)
        k in keys(fm.design) || throw(ArgumentError(
            "unknown design variable :$k (declared: $(join(keys(fm.design), ", ")))"))
        v isa Real || throw(ArgumentError("design values must be scalar Reals"))
        fm.design = merge(fm.design, NamedTuple{(k,)}((Float64(v),)))
    end
    for (k, θ) in pairs(fm.design)
        _write_design!(fm.slots[k], θ)
    end
    fm.B = nothing; fm.J = nothing
    return fm.design
end

"""
    value_and_gradient(fm::FrequencyMatch) -> (J, G)

[`run_forward!`](@ref) followed by [`gradient!`](@ref) — the single protocol
pair optimizers drive.
"""
function value_and_gradient(fm::FrequencyMatch)
    J = run_forward!(fm)
    return J, gradient!(fm)
end
