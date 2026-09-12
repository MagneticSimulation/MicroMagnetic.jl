# ---------------------------------------------------------------------------
# SweepMatch — stage-chain inverse-design scenario.
#
# Forward: a per-stage Zeeman protocol; every stage is one steady relax that
# continues from the previous stage's terminal state (path continuity, the
# sim_with stage-chain semantics).  Each stage's converged state m*_k and its
# total field H0_k are checkpointed (~len(stages) × 6N Float64) together with
# the observed curve.
#
# Backward: the steady-state adjoint per stage, traversed in REVERSE
# stage order.  For a stage loss g_k the adjoint field λ_k solves
#
#     (Df)ᵀ λ_k = P_t ∇_m g_k(m*_k)        (Df from linearize.jl at (m*_k, H0_k))
#
# and contributes  ∂J/∂θ += r_kᵀ ∂H/∂θ  with the theta weight
# r_k = γ′(λ_k×m*_k − α m*_k×(λ_k×m*_k))  (llg_param_weight_kernel!), applied
# to the design variables through the grad_* kernels.  The λ "chain"
# between stages is empty by construction: every stage is a steady state, so
# stage k's gradient only involves stage k's checkpoint — the path dependence
# lives entirely in the checkpointed states.  A stage with zero loss weight
# has P_t∇g = 0 and solve_steady_adjoint returns λ = 0 immediately.
#
# Why the LLG linearization is valid for SD-driven stages: every driver in the
# package relaxes onto the same equilibrium manifold {m×H_eff = 0} on which
# f_LL = −γ′(m×H + α m×(m×H)) also vanishes, so the implicit derivative
# dm*/dθ = −(Df)⁻¹ ∂f/∂θ — and hence the adjoint gradient — is independent of
# which flow reached the state.  α only enters the adjoint solve's
# conditioning (α > 0 makes the tangent operator dissipative); the gradient
# is invariant to it, and to the forward driver's parameters.
#
# Scope guards (all fail fast): MicroSim on FDMesh, Float64 state, full
# occupancy (no vacuum cells — the steady solver requires unit m per site),
# no pinned spins, design variables :Ku/:Aex (the grad kernels wired so far),
# interactions limited to those with a Kᵀ in linearize.jl (throw list there).
# ---------------------------------------------------------------------------

export Moment, CurveTarget, SquaredError, SweepMatch, run_forward!, gradient!,
       set_design!, value_and_gradient, checkpoint_memory

"""
    Moment(comp::Symbol)

Observation handle for `CurveTarget`: the Ms-weighted average of one
magnetization component at a stage,

    M_comp = Σ_i mu0_Ms_i · m_i,comp / Σ_i mu0_Ms_i ,

i.e. the magnetic moment per unit volume (the μ0 factor cancels).  Vacuum
cells drop out through their zero `mu0_Ms`.  `comp` is `:x`, `:y` or `:z`.
"""
struct Moment
    comp::Symbol
    function Moment(comp::Symbol)
        comp in (:x, :y, :z) ||
            throw(ArgumentError("Moment expects :x, :y or :z (got $comp)"))
        return new(comp)
    end
end

_comp_index(c::Symbol) = c === :x ? 1 : (c === :y ? 2 : 3)

"""
    CurveTarget(obs::Moment, values; weights=nothing)

Target curve for `SweepMatch`: `values` holds one target value per stage (a
single value is broadcast to all stages); `weights` holds the per-stage loss
weights (default 1 per stage, a single value broadcasts).
"""
struct CurveTarget
    obs::Moment
    values::Vector{Float64}
    weights::Union{Nothing,Vector{Float64}}
end
CurveTarget(obs::Moment, values; weights = nothing) = CurveTarget(obs, values, weights)

# One protocol stage: a static Zeeman field in A/m.  A plain struct so P4 can
# widen it (dynamic stages, per-stage driver) without changing the sweep API.
struct _SweepStage
    H::NTuple{3,Float64}
end

# A declared design variable: which parameter, on which interaction, and the
# gradient shape the user declared (scalar ↔ Fill, map ↔ dense array).
struct _DesignVar
    name::Symbol
    kind::Symbol            # :Ku | :Aex
    inter::Any              # Anisotropy / Exchange inside sim.interactions
    is_scalar::Bool
end

mutable struct SweepMatch
    sim::AbstractSim
    stages::Vector{_SweepStage}
    target::CurveTarget
    loss::SquaredError
    design::NamedTuple
    G_maps::NamedTuple            # per-cell accumulators, one per design var
    zeeman_name::String

    # forward options
    driver::String
    stopping_dmdt::Float64
    max_steps::Int
    relax_tol::Float64
    relax_residual::Vector{Float64}   # max |m×H| (A/m) per stage after relax

    # adjoint options (α/γ are the *linearization* parameters — see module head)
    alpha::Float64
    gamma::Float64
    adjoint_backend::Symbol
    adjoint_tol::Float64
    adjoint_maxiter::Int

    # state
    forward_ready::Bool
    J::Float64
    G::Any                        # NamedTuple with the declared gradient shapes
    m0::Vector{Float64}           # initial state recorded at construction
    obs_curve::Vector{Float64}    # observed value per stage
    cp_m::Vector{Vector{Float64}} # checkpointed states m*_k (unit vectors)
    cp_H0::Vector{Vector{Float64}}# checkpointed total fields H0_k at m*_k
    adjoint_residual::Vector{Float64}
    adjoint_bnorm::Vector{Float64}  # ‖P_t ∇g‖ per stage (degenerate ⇒ ~0)

    # scratch
    mu0_Ms_host::Vector{Float64}
    ms_sum::Float64
    obs_comp::Int
    ws::AdjointWorkspace
    grad_g::Vector{Float64}
    r::Vector{Float64}
end

# ---- sweep protocol parsing -------------------------------------------------

function _parse_sweep(sweep::NamedTuple)
    iset = keys(sweep)
    iset ⊆ (:H_list, :angles, :directions) ||
        throw(ArgumentError("sweep accepts the keys :H_list with one of :angles / " *
                            ":directions (got $iset)"))
    haskey(sweep, :H_list) ||
        throw(ArgumentError("sweep needs :H_list (field magnitude(s) in A/m)"))
    haskey(sweep, :angles) && haskey(sweep, :directions) &&
        throw(ArgumentError("sweep accepts either :angles or :directions, not both"))

    if haskey(sweep, :angles)
        angles = collect(Float64, sweep.angles)
        n = length(angles)
        n > 0 || throw(ArgumentError("sweep angles must be non-empty"))
        Hs = _expand_H(sweep.H_list, n)
        # documented convention: rotate in the x–z plane through the easy axis,
        # Ĥ(φ) = (sin φ, 0, cos φ); use :directions for any other protocol
        return [_SweepStage((H * sin(φ), 0.0, H * cos(φ))) for (H, φ) in zip(Hs, angles)]
    else
        haskey(sweep, :directions) ||
            throw(ArgumentError("sweep needs :angles or :directions alongside :H_list"))
        dirs = sweep.directions
        n = length(dirs)
        n > 0 || throw(ArgumentError("sweep directions must be non-empty"))
        Hs = _expand_H(sweep.H_list, n)
        stages = _SweepStage[]
        for (H, d) in zip(Hs, dirs)
            length(d) == 3 || throw(ArgumentError("each direction must have 3 components"))
            nrm = sqrt(Float64(d[1])^2 + Float64(d[2])^2 + Float64(d[3])^2)
            nrm > 0 || throw(ArgumentError("sweep directions must be nonzero"))
            push!(stages, _SweepStage((H * d[1] / nrm, H * d[2] / nrm, H * d[3] / nrm)))
        end
        return stages
    end
end

function _expand_H(H_list, n::Int)
    Hs = H_list isa Number ? fill(Float64(H_list), n) : Float64.(collect(H_list))
    length(Hs) == n ||
        throw(DimensionMismatch("H_list has length $(length(Hs)) but the sweep has " *
                                "$n stages"))
    return Hs
end

# ---- design parsing ---------------------------------------------------------

function _parse_design(sim::AbstractSim, design::NamedTuple)
    isempty(design) && throw(ArgumentError("design must declare at least one variable, " *
                                           "e.g. design = (Ku = Ku_map)"))
    vars = []
    for (k, v) in pairs(design)
        k in (:Ku, :Aex) ||
            throw(ArgumentError("design variable :$k is not wired yet (supported: " *
                                ":Ku, :Aex; :Ms/:D follow with their grad kernels)"))
        if k === :Ku
            inter = findfirst(i -> i isa Anisotropy, sim.interactions)
            inter === nothing &&
                throw(ArgumentError("design :Ku needs a uniaxial anisotropy: call " *
                                    "add_anis(sim, Ku; axis=...) before SweepMatch"))
            inter = sim.interactions[inter]
            kind = :Ku
        else
            inter = findfirst(i -> i isa Exchange, sim.interactions)
            inter === nothing &&
                throw(ArgumentError("design :Aex needs an exchange interaction: call " *
                                    "add_exch(sim, A) before SweepMatch"))
            inter = sim.interactions[inter]
            kind = :Aex
        end
        arr = make_param(Float64, v, sim.mesh, sim.n_total)
        push!(vars, _DesignVar(k, kind, inter, arr isa Fill))
    end
    return NamedTuple{keys(design)}(Tuple(vars))
end

function _apply_design!(prob::SweepMatch, var::_DesignVar, v)
    sim = prob.sim
    if var.kind === :Ku
        # replace (a Fill cannot be written in place); anisotropy has no caches
        var.inter.Ku = make_param(Float64, v, sim.mesh, sim.n_total)
    else
        arr = make_param(Float64, v, sim.mesh, sim.n_total)
        var.inter.Ax = arr
        var.inter.Ay = arr
        var.inter.Az = arr
        # the stencil pair tables cache harmonic means of the OLD stiffness —
        # dropping them forces a rebuild from the new arrays
        var.inter.pair_Ax = var.inter.pair_Ay = var.inter.pair_Az = nothing
        var.inter.stencil_layout = -1
    end
    return var
end

# ---- constructor ------------------------------------------------------------

"""
    SweepMatch(sim; design, sweep, target, loss=SquaredError(), driver="SD",
               stopping_dmdt=1e-9, max_steps=20000, relax_tol=1e-4,
               zeeman_name="sweep_H", alpha=0.05, gamma=2.21e5,
               adjoint_backend=:krylov, adjoint_tol=1e-8, adjoint_maxiter=10_000)

Stage-chain inverse problem: match an
observed curve — e.g. a rotation magnetization curve — by sweeping a Zeeman
protocol over steady stages and comparing against `target`.

- `design`: NamedTuple of design variables, e.g. `design = (Ku = Ku_map)` or
  `design = (Aex = 1.3e-11, Ku = Ku_map)`.  A `Number` declares a scalar
  variable (gradient returned as a scalar); an array/function declares a
  spatial map (gradient returned as a per-cell vector).  Supported: `:Ku`
  (uniaxial anisotropy, needs `add_anis` first) and `:Aex` (exchange, needs
  `add_exch` first).  The declared values are pushed onto `sim` at construction.
- `sweep`: the stage protocol, one field per stage (A/m):
  `(H_list = 80e3, angles = collect(range(0, 2π; length=25)))` rotates the
  field in the **x–z plane**, `Ĥ(φ) = (sin φ, 0, cos φ)` — through the easy
  axis, which is what makes a `Moment` curve respond to `Ku`.  For any other
  protocol pass `(H_list = ..., directions = [(hx, hy, hz), ...])` (directions
  are normalized).  `H_list` may be a number (constant magnitude) or one value
  per stage.
- `target`: `CurveTarget(Moment(:z), values; weights=...)` with one target
  value per stage; every stage with a nonzero weight contributes its own
  steady adjoint (zero-weight stages contribute nothing — their λ is the
  projection of the zero vector, i.e. 0).
- `loss`: `SquaredError()` (only variant so far).
- `driver`: forward relaxation driver per stage, `"SD"` (default) or `"LLG"`.
  SD is strongly recommended for stage chains: at equal `stopping_dmdt` it
  lands orders of magnitude closer to the equilibrium (its stopping measure is
  step length, not precession rate).
- `stopping_dmdt`, `max_steps`: forwarded to `relax` at every stage.  The
  default `1e-9` is adjoint-grade: the gradient inherits the relaxation
  residual, so tight stages are what makes `gradient!` match finite
  differences.  `relax_tol` (default `1e-4` A/m) is the post-relax quality
  gate: SD stages whose max |m×H| exceeds it are re-relaxed with a fresh
  Barzilai–Borwein phase (up to 3 attempts) because a converged previous
  stage can collapse the BB step under the new field; per-stage residuals
  land in `prob.relax_residual`.
- `alpha`, `gamma`: damping and gyromagnetic ratio of the **adjoint
  linearization** (not the forward driver's).  The gradient is invariant to
  both (module head comment); α > 0 conditions the adjoint solve.
- `adjoint_backend`: `:krylov` (GMRES, default) or `:pseudotime` — the two
  `solve_steady_adjoint` backends.

The sim must be Float64, fully occupied (no vacuum), unpinned, initialized
with `init_m0`, and must not already carry a Zeeman named `zeeman_name`.
`run_forward!` always restarts from the state recorded at construction, so
`J(θ)` is a deterministic function of the design.

Protocol (`run_forward!` → `gradient!` → `set_design!` → …): `set_design!`
invalidates the forward; `gradient!` must follow `run_forward!` at the same
design values.
"""
function SweepMatch(sim::AbstractSim; design, sweep, target, loss = SquaredError(),
                    driver::AbstractString = "SD",
                    stopping_dmdt::Real = 1e-9,
                    max_steps::Integer = 20000,
                    relax_tol::Real = 1e-4,
                    zeeman_name::AbstractString = "sweep_H",
                    alpha::Real = 0.05,
                    gamma::Real = 2.21e5,
                    adjoint_backend::Symbol = :krylov,
                    adjoint_tol::Real = 1e-8,
                    adjoint_maxiter::Integer = 10_000)
    sim isa MicroSim ||
        throw(ArgumentError("SweepMatch supports MicroSim (FDMesh) only"))
    eltype(sim.spin) === Float64 ||
        throw(ArgumentError("the adjoint layer is Float64-only (sim.spin is " *
                            "$(eltype(sim.spin))); build the sim with default precision"))
    any(Array(sim.spin) .!= 0) ||
        throw(ArgumentError("sim.spin is zero: call init_m0(sim, m0) before SweepMatch"))
    any(sim.pins) &&
        throw(ArgumentError("pinned spins are not supported by the adjoint layer"))
    driver in ("SD", "LLG") ||
        throw(ArgumentError("driver must be \"SD\" or \"LLG\" (got $driver)"))
    adjoint_backend in (:krylov, :pseudotime) ||
        throw(ArgumentError("adjoint_backend must be :krylov or :pseudotime"))

    loss isa SquaredError ||
        throw(ArgumentError("only loss = SquaredError() is supported so far"))

    stages = _parse_sweep(sweep)
    n = length(stages)

    # full occupancy: the steady adjoint solver requires unit vectors per site
    ms_host = Float64.(Array(sim.mu0_Ms))
    all(ms_host .> 0) ||
        throw(ArgumentError("SweepMatch requires a fully occupied mesh (vacuum cells, " *
                            "Ms = 0, are not supported by the steady adjoint solver yet)"))

    # normalize target lengths (1 broadcasts to n)
    tvals = Float64.(collect(target.values))
    (length(tvals) == n || length(tvals) == 1) ||
        throw(DimensionMismatch("target.values must have one entry per stage ($n) " *
                                "or a single entry, got $(length(tvals))"))
    length(tvals) == 1 && (tvals = fill(tvals[1], n))
    tw = target.weights === nothing ? ones(n) : Float64.(collect(target.weights))
    (length(tw) == n || length(tw) == 1) ||
        throw(DimensionMismatch("target.weights must have one entry per stage ($n) " *
                                "or a single entry"))
    length(tw) == 1 && (tw = fill(tw[1], n))
    target = CurveTarget(target.obs, tvals, tw)

    dvars = _parse_design(sim, design)
    G_maps = NamedTuple{keys(dvars)}(Tuple(zeros(Float64, sim.n_total) for _ in pairs(dvars)))

    # the sweep owns one dedicated Zeeman interaction
    any(inter -> inter.name == zeeman_name, sim.interactions) &&
        throw(ArgumentError("sim already has an interaction named \"$(zeeman_name)\"; " *
                            "pass a different zeeman_name"))

    prob = SweepMatch(sim, stages, target, loss, dvars, G_maps, String(zeeman_name),
                      String(driver), Float64(stopping_dmdt), Int(max_steps),
                      Float64(relax_tol), zeros(n),
                      Float64(alpha), Float64(gamma), adjoint_backend,
                      Float64(adjoint_tol), Int(adjoint_maxiter),
                      false, NaN, nothing,
                      Float64.(Array(sim.spin)), zeros(n), Vector{Vector{Float64}}(undef, n),
                      Vector{Vector{Float64}}(undef, n), zeros(n), zeros(n),
                      ms_host, sum(ms_host), _comp_index(target.obs.comp),
                      AdjointWorkspace(sim), zeros(3 * sim.n_total), zeros(3 * sim.n_total))

    # push the declared design onto the sim and prepare driver + protocol field
    for (k, v) in pairs(design)
        _apply_design!(prob, dvars[k], v)
    end
    set_driver(sim; driver = prob.driver)
    add_zeeman(sim, stages[1].H; name = prob.zeeman_name)
    return prob
end

# ---- helpers ----------------------------------------------------------------

function _normalize_sites!(m::Vector{Float64})
    @inbounds for i in 1:3:length(m)
        n = sqrt(m[i]^2 + m[i + 1]^2 + m[i + 2]^2)
        if n > 0 && n != 1.0
            m[i] /= n
            m[i + 1] /= n
            m[i + 2] /= n
        end
    end
    return m
end

function _weighted_moment(m::Vector{Float64}, ms::Vector{Float64}, comp::Int,
                          ms_sum::Float64)
    s = 0.0
    @inbounds for i in 1:length(ms)
        s += ms[i] * m[3 * (i - 1) + comp]
    end
    return s / ms_sum
end

"""
    checkpoint_memory(prob::SweepMatch) -> Int

Bytes held by the stage checkpoints (`m*_k` + `H0_k`, ~stages × 6N Float64;
the observed curve and residuals are negligible and not counted).
"""
checkpoint_memory(prob::SweepMatch) =
    sum(sizeof, prob.cp_m; init = 0) + sum(sizeof, prob.cp_H0; init = 0)

# max |m×H| over sites (A/m) — host-side relaxation-quality diagnostic.
function _torque_residual(sim::AbstractSim)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    m = Array(sim.spin)
    h = Array(sim.field)
    tmax = 0.0
    @inbounds for i in 1:sim.n_total
        j = 3 * (i - 1)
        tx = m[j + 2] * h[j + 3] - m[j + 3] * h[j + 2]
        ty = m[j + 3] * h[j + 1] - m[j + 1] * h[j + 3]
        tz = m[j + 1] * h[j + 2] - m[j + 2] * h[j + 1]
        tmax = max(tmax, sqrt(tx^2 + ty^2 + tz^2))
    end
    return tmax
end

# One steady stage for the SD driver: a fresh Barzilai–Borwein phase per stage
# (the BB rule degenerates when seeded with the ~zero displacement of a
# converged previous stage under a new field — the step collapses and relax
# stops without moving), plus a torque-residual guard with bounded retries.
function _relax_stage!(prob::SweepMatch)
    sim = prob.sim
    if sim.driver isa SD
        res = Inf
        for _ in 0:2
            sim.driver.steps = 0
            relax(sim; max_steps = prob.max_steps, stopping_dmdt = prob.stopping_dmdt,
                  save_data_every = -1, save_m_every = -1)
            res = _torque_residual(sim)
            res <= prob.relax_tol && return res
        end
        @warn "SweepMatch stage did not reach relax_tol after 3 relax attempts " *
              "(max |m×H| = $(res) A/m); the gradient for this stage inherits the error"
        return res
    else
        relax(sim; max_steps = prob.max_steps, stopping_dmdt = prob.stopping_dmdt,
              save_data_every = -1, save_m_every = -1)
        return _torque_residual(sim)
    end
end

# ---- forward ----------------------------------------------------------------

"""
    run_forward!(prob::SweepMatch) -> Float64

Run the stage chain: restart from the construction-time state, then for each
stage update the sweep field, relax to steady state (continuing from the
previous stage's terminal state), and checkpoint `m*_k`, `H0_k` and the
observed value.  SD stages restart the Barzilai–Borwein phase and are verified
against `relax_tol` (max |m×H|, with bounded retries — see `_relax_stage!`).
Returns `J = Σ_k w_k (M_k − target_k)²`.
"""
function run_forward!(prob::SweepMatch)
    sim = prob.sim
    copyto!(sim.spin, prob.m0)                       # deterministic restart
    hasproperty(sim.driver, :integrator) &&
        set_initial_condition!(sim, sim.driver.integrator)

    J = 0.0
    for k in eachindex(prob.stages)
        update_zeeman(sim, prob.stages[k].H; name = prob.zeeman_name)
        prob.relax_residual[k] = _relax_stage!(prob)
        m = Float64.(Array(sim.spin))
        _normalize_sites!(m)
        copyto!(sim.spin, m)                          # checkpoint ⇄ state consistency
        prob.cp_m[k] = m
        MicroMagnetic.effective_field(sim, sim.spin, 0.0)
        prob.cp_H0[k] = Float64.(Array(sim.field))
        Mk = _weighted_moment(m, prob.mu0_Ms_host, prob.obs_comp, prob.ms_sum)
        prob.obs_curve[k] = Mk
        J += prob.target.weights[k] * (Mk - prob.target.values[k])^2
    end
    prob.J = J
    prob.forward_ready = true
    return J
end

# ---- backward ---------------------------------------------------------------

# Projected relative residual ‖P_t(Df)ᵀλ − P_t∇g‖ / ‖P_t∇g‖ — the quantity the
# solvers actually target (the unprojected (Df)ᵀλ has a benign longitudinal
# part: Df·m̂ ≠ 0 once anisotropy/demag break rotational equivariance, so the
# raw residual is dominated by that leakage, not by the solve).
# Returns (relres, ‖P_t∇g‖); relres = 0.0 when the right-hand side is below the
# solver's zero short-circuit (λ = 0 exactly, no solve).  For a stage whose
# ‖P_t∇g‖ sits barely above that threshold the RELATIVE residual amplifies
# roundoff and is meaningless — judge such stages by their absolute size
# (‖P_t∇g‖ itself bounds the error contribution) via the returned norm.
function _stage_residual(vjp!::Function, lam::Vector{Float64}, grad_g::Vector{Float64},
                         m_star::Vector{Float64})
    b = project_tangent(grad_g, m_star)
    nb = norm(b)
    nb <= 128 * eps(Float64) * max(norm(grad_g), 1.0) && return 0.0, nb
    out = zeros(length(lam))
    vjp!(out, lam)
    project_tangent!(out, out, m_star)
    return norm(out .- b) / nb, nb
end

"""
    gradient!(prob::SweepMatch) -> NamedTuple

Reverse-stage-order steady adjoint: for k = nstages,…,1 solve
`(Df)ᵀλ_k = P_t∇g_k` at the checkpointed `(m*_k, H0_k)`, form the theta weight
`r_k = γ′(λ_k×m*_k − α m*_k×(λ_k×m*_k))`, and accumulate `∂J/∂θ += r_kᵀ ∂H/∂θ`
with the grad kernels.  Returns a NamedTuple mirroring `design`: scalar
variables collapse the per-cell gradient by summation, map variables return
the per-cell vector (a fresh copy each call).  Requires a fresh
`run_forward!` (a `set_design!` in between invalidates it).

Per-stage adjoint diagnostics are stored in `prob.adjoint_residual` (relative
projected residual; `0.0` marks stages whose loss weight is zero — their λ is
exactly zero without any solve) and `prob.adjoint_bnorm` (‖P_t∇g‖; a stage with
‖P_t∇g‖ ≈ 0 contributes nothing to the gradient regardless of its residual).
"""
function gradient!(prob::SweepMatch)
    prob.forward_ready ||
        throw(ErrorException("gradient! needs a fresh forward: call run_forward!(prob) " *
                             "after the last set_design!/construction"))
    sim = prob.sim
    n_total = sim.n_total
    back = KernelAbstractions.CPU()
    ap = Float64(prob.alpha)
    gp = Float64(prob.gamma)
    gp_eff = gp / (1 + ap^2)
    comp = prob.obs_comp

    for G in prob.G_maps
        fill!(G, 0.0)
    end

    for k in length(prob.stages):-1:1                # strict reverse stage order
        m_star = prob.cp_m[k]
        H0 = prob.cp_H0[k]

        # ∇_m g_k = 2 w_k (M_k − c_k) · (mu0_Ms_i / Σ mu0_Ms) ê_comp
        fill!(prob.grad_g, 0.0)
        gk = 2 * prob.target.weights[k] * (prob.obs_curve[k] - prob.target.values[k]) /
             prob.ms_sum
        if gk != 0.0
            @inbounds for i in 1:n_total
                prob.grad_g[3 * (i - 1) + comp] = gk * prob.mu0_Ms_host[i]
            end
        end

        vjp! = function (out, x)
            return llg_rhs_vjp!(out, sim, m_star, H0, x, prob.ws; alpha = ap, gamma = gp)
        end
        lam = solve_steady_adjoint(vjp!, prob.grad_g, m_star;
                                   backend = prob.adjoint_backend,
                                   tol = prob.adjoint_tol,
                                   maxiter = prob.adjoint_maxiter)
        prob.adjoint_residual[k], prob.adjoint_bnorm[k] =
            _stage_residual(vjp!, lam, prob.grad_g, m_star)

        # theta weight r = γ′(λ×m − α m×(λ×m));  ∂J/∂θ = rᵀ ∂H/∂θ
        llg_param_weight_kernel!(back, groupsize[])(prob.r, m_star, lam, ap, gp_eff;
                                                    ndrange = n_total)

        for name in keys(prob.design)
            var = prob.design[name]
            G = prob.G_maps[name]
            if var.kind === :Ku
                grad_ku_kernel!(back, groupsize[])(G, m_star, prob.r,
                                                   var.inter.axis_x, var.inter.axis_y,
                                                   var.inter.axis_z, sim.mu0_Ms;
                                                   ndrange = n_total)
            else
                mesh = sim.mesh
                grad_aex_kernel!(back, _stencil_wg(mesh))(G, m_star, prob.r,
                    sim.mu0_Ms, sim.inv_ms, var.inter.Ax, var.inter.Ay, var.inter.Az,
                    Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                    mesh.nx, mesh.ny, mesh.nz,
                    mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                    ndrange = (mesh.nx, mesh.ny, mesh.nz))
            end
        end
    end

    prob.G = NamedTuple{keys(prob.design)}(Tuple(var.is_scalar ? sum(prob.G_maps[name]) :
                                                 copy(prob.G_maps[name])
                                                 for (name, var) in pairs(prob.design)))
    return prob.G
end

# ---- iteration protocol -----------------------------------------------------

"""
    set_design!(prob::SweepMatch; kwargs...)

Write new design values back into the sim (Ku → the `Anisotropy` interaction,
Aex → the `Exchange` interaction, dropping its stencil pair tables) and
invalidate the forward: the next `gradient!` requires a fresh `run_forward!`.
Each keyword must be a declared design variable; scalars for scalar
variables, arrays/functions for map variables.
"""
function set_design!(prob::SweepMatch; kwargs...)
    isempty(kwargs) && throw(ArgumentError("set_design! needs at least one design " *
                                           "variable, e.g. set_design!(prob; Ku = map)"))
    for (k, v) in kwargs
        haskey(prob.design, k) ||
            throw(ArgumentError("unknown design variable :$k (declared: " *
                                "$(keys(prob.design)))"))
        var = prob.design[k]
        if var.is_scalar && !(v isa Number)
            throw(ArgumentError("design :$k was declared as a scalar; set_design! " *
                                "expects a Number"))
        elseif !var.is_scalar && v isa Number
            throw(ArgumentError("design :$k was declared as a spatial map; set_design! " *
                                "expects an array or function"))
        end
        _apply_design!(prob, var, v)
    end
    prob.forward_ready = false
    return prob
end

"""
    value_and_gradient(prob::SweepMatch) -> (J, G)

`run_forward!` followed by `gradient!` — the single entry point optimizers
need.
"""
function value_and_gradient(prob::SweepMatch)
    return run_forward!(prob), gradient!(prob)
end
