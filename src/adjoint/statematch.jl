# ---------------------------------------------------------------------------
# StateMatch: static-configuration inverse design.
#
#   sm = StateMatch(sim; design=(Ku=..., Aex=..., D=...), target_m=m_target)
#   J = run_forward!(sm)      # write design → fixed init_m0 → relax → freeze m*, H0
#   G = gradient!(sm)         # one steady-state adjoint solve + ∂f/∂θ accumulation
#   set_design!(sm; ...)      # write new design values back (with cache drops)
#
# Forward: one relaxation to a steady state m* (f(m*) = 0) from a fixed initial
# state — no sweep, no path dependence.  Backward: the single-stage
# steady-state adjoint, (Df)ᵀλ = P_t ∇g, solved by steady.jl against the
# llg_rhs_vjp! VJP; then ∂J/∂θ = w·(∂H/∂θ) with the reduced weight
# w = γ′(λ×m* − α m*×(λ×m*)) (sign and γ′ folded in) accumulated by the
# grad_ku/grad_aex/grad_dmi kernels.  The FD tests (test_statematch_fd.jl) are
# the final judge of every sign/factor convention.
#
# Float64/CPU only at this stage (the whole adjoint layer is Vector{Float64});
# the design write-back reuses the existing parameter machinery (make_param,
# set_Ms) so Fill designs stay O(1) and gradients auto-merge to scalars.
# ---------------------------------------------------------------------------

export StateMatch, SquaredError, run_forward!, gradient!, set_design!,
       value_and_gradient

"""
    SquaredError()
    SquaredError(weight)

Per-cell squared-error metric for `StateMatch`: the loss is

    J(m*) = Σ_i  w_i · |m*_i − target_i|²

with `w_i = 1` for the empty constructor.  `weight` is a length-`n_total`
vector of nonnegative per-cell weights (e.g. fit only a sub-region by setting
zeros elsewhere).  Vacuum cells (`Ms == 0`) are always excluded from the loss.
"""
struct SquaredError
    weight::Vector{Float64}     # per-cell; empty = uniform weight 1
end

SquaredError() = SquaredError(Float64[])

function SquaredError(weight::AbstractVector{<:Real})
    all(w -> w >= 0, weight) ||
        throw(ArgumentError("SquaredError weights must be nonnegative"))
    return SquaredError(collect(Float64, weight))
end

# Design keys supported by this scenario (θ enters f only through H; the Ms
# gradient needs the demag linear pass-through and is not yet implemented).
const _STATEMATCH_KEYS = (:Ku, :Aex, :D)

mutable struct StateMatch
    sim::AbstractSim
    design::NamedTuple                    # declared design values (Number/Fill/array/function)
    target_m::Vector{Float64}             # 3N unit-vector target (vacuum cells zero)
    metric::SquaredError
    init_m::TupleOrArrayOrFunction        # fixed initial state for every forward
    anis::Union{Nothing,Anisotropy}
    exch::Union{Nothing,Exchange}
    dmi::Union{Nothing,DMI}
    # frozen forward state (filled by run_forward!)
    m_star::Union{Nothing,Vector{Float64}}
    H0::Union{Nothing,Vector{Float64}}
    ws::Union{Nothing,AdjointWorkspace}
    J::Float64
    max_transverse::Float64               # max |m*×H0| — convergence diagnostic
    grad::Union{Nothing,NamedTuple}
    # options
    max_steps::Int
    stopping_dmdt::Float64
    adjoint_backend::Symbol
    adjoint_tol::Float64
    adjoint_maxiter::Int
    verbose::Bool
end

"""
    StateMatch(sim; design, target_m, metric=SquaredError(), init_m=(0, 0, 1),
               max_steps=10000, stopping_dmdt=1e-8, backend=:krylov,
               tol=1e-10, maxiter=10_000, verbose=false)

Static-configuration inverse design: find
the material maps `θ` that drive the relaxed state `m*(θ)` towards a target
configuration.

- `design`: which parameters are design variables, e.g. `design = (Ku = 5e4,
  Aex = 1.3e-11, D = 3e-3)`.  Each value may be a scalar (uniform parameter),
  an `n_total`-length array, a spatial function `f(i,j,k,dx,dy,dz)` / `f(x,y,z)`
  / `f(region_id)` (e.g. [`region_map`](@ref)), or an `O(1)` `Fill`.  The
  gradient returned by [`gradient!`](@ref) matches the declared shape: scalars
  and `Fill`s give a scalar, arrays/functions give the per-cell map.
- `target_m`: the 3N target configuration — an array, a tuple or a spatial
  function (same forms as [`init_m0`](@ref)).  It is normalised per site, and
  vacuum cells (`Ms == 0`) are zeroed and excluded from the loss.
- `metric`: [`SquaredError`](@ref), optionally with per-cell weights.
- `init_m`: the fixed initial state every forward relaxes from (same forms as
  `init_m0`; a deterministic inverse problem needs a fixed start).
- `max_steps`, `stopping_dmdt`: forwarded to [`relax`](@ref).  The gradient is
  only as good as the forward convergence — a warning is printed when the
  relaxation stops too loose for trustworthy gradients (see
  `sm.max_transverse`).
- `backend`, `tol`, `maxiter`: the steady-state adjoint solver
  ([`solve_steady_adjoint`](@ref)) options.  The default `:krylov` (GMRES on
  `P_t (Df)ᵀ P_t`) is the robust choice; `:pseudotime` additionally requires
  the symmetric part of `(Df)ᵀ` to be negative definite on the tangent space,
  which fails at stiffness-compensated equilibria (e.g. a Zeeman field
  balancing the anisotropy field — one tangential curvature nearly vanishes
  and the precession coupling makes the symmetric part indefinite).

Forward-convergence guidance: with the LLG driver the fast precession about
the (non-vanishing) longitudinal field pins the adaptive step near equilibrium,
so the `dmdt` tail converges very slowly; `sim.driver.precession = false`
(pure damping — same equilibrium branch, same gradient) removes the pinning.
The `dmdt` stopping floor also sits at the integrator's noise plateau: for
`stopping_dmdt ≤ 1e-8`, tighten the RK tolerance as well, e.g.
`sim.driver.integrator.tol = 1e-13` (the default `1e-6` stalls the residual
around `|m×H| ~ 50 A/m`).

The forward must run with the `LLG` driver and a uniform damping α (the
adjoint linearises exactly this flow and reuses the driver's α and γ); heavy
damping (e.g. `set_alpha(sim, 1.0)`) makes the relaxations much cheaper.

Requires `add_anis` / `add_exch` / `add_dmi` interactions for every declared
design key (the first uniaxial anisotropy / exchange / DMI in `sim.interactions`
is used).  Supported interactions: Exchange, uniaxial Anisotropy, bulk and
interfacial DMI, Zeeman, Demag — everything else (stochastic fields, torques,
...) is rejected by the VJP layer.
"""
function StateMatch(sim::AbstractSim; design::NamedTuple = NamedTuple(),
                    target_m, metric::SquaredError = SquaredError(),
                    init_m::TupleOrArrayOrFunction = (0, 0, 1),
                    max_steps::Int = 10000, stopping_dmdt::Float64 = 1e-8,
                    backend::Symbol = :krylov, tol::Float64 = 1e-10,
                    maxiter::Int = 10_000, verbose::Bool = false)
    sim.driver isa LLG || throw(ArgumentError(
        "StateMatch relaxes with the LLG driver (the adjoint linearises this " *
        "flow); got $(nameof(typeof(sim.driver))). Use `set_driver(sim; driver=\"LLG\")` " *
        "and `set_alpha` first."))
    _adjoint_alpha_gamma(sim)   # fail fast on spatial α

    isempty(design) && throw(ArgumentError(
        "StateMatch needs at least one design key among $(_STATEMATCH_KEYS)"))
    for k in keys(design)
        k in _STATEMATCH_KEYS || throw(ArgumentError(
            "unsupported design key :$k (supported: $(_STATEMATCH_KEYS); " *
            "Ms design needs the not-yet-implemented ∂H/∂Ms machinery)")
        )
    end

    N = sim.n_total
    anis = :Ku in keys(design) ? _find_interaction(Anisotropy, sim, :Ku) : nothing
    exch = :Aex in keys(design) ? _find_interaction(Exchange, sim, :Aex) : nothing
    dmi = :D in keys(design) ? _find_interaction(DMI, sim, :D) : nothing

    sm = StateMatch(sim, design, _materialise_target(sim, target_m), metric, init_m,
                    anis, exch, dmi, nothing, nothing, nothing, NaN, NaN, nothing,
                    max_steps, stopping_dmdt, backend, tol, maxiter, verbose)
    _apply_design!(sm)
    return sm
end

# first interaction of the given type; nothing found → helpful error
function _find_interaction(::Type{I}, sim::AbstractSim, key::Symbol) where {I<:MicroEnergy}
    for inter in sim.interactions
        inter isa I && return inter
    end
    return throw(ArgumentError(
        "design key :$key needs a $(nameof(I)) interaction; add one " *
        "(e.g. $(nameof(I) === :Anisotropy ? "add_anis" : nameof(I) === :Exchange ? "add_exch" : "add_dmi")(sim, ...)) first"))
end

# The adjoint VJP reuses the driver's α and γ; spatial α would need a kernel
# variant, so require the uniform (Fill) representation here.
function _adjoint_alpha_gamma(sim::AbstractSim)
    driver = sim.driver
    driver isa LLG || throw(ArgumentError("StateMatch requires the LLG driver"))
    driver.alpha isa Fill || throw(ArgumentError(
        "StateMatch requires uniform damping (alpha stored as a Fill, set via " *
        "`set_alpha(sim, <number>)`); spatial alpha is not supported by the " *
        "pointwise VJP kernel yet."))
    return Float64(driver.alpha.value), Float64(driver.gamma)
end

# Materialise the target spec: 3N Float64, normalised per site, zero at vacuum.
function _materialise_target(sim::AbstractSim, target)
    N = sim.n_total
    t = zeros(Float64, 3N)
    init_vector!(t, sim.mesh, target)
    normalise(t, N)
    ms = Array(sim.mu0_Ms)
    @inbounds for i in 1:N
        ms[i] == 0 && (t[3i-2] = t[3i-1] = t[3i] = 0.0)
    end
    return t
end

# per-cell loss weights: metric weight × vacuum mask (no factor 2 here)
function _loss_weights(sm::StateMatch)
    N = sm.sim.n_total
    w = zeros(Float64, N)
    weight = sm.metric.weight
    ms = Array(sm.sim.mu0_Ms)
    @inbounds for i in 1:N
        wi = isempty(weight) ? 1.0 : weight[i]
        w[i] = ms[i] == 0 ? 0.0 : wi
    end
    return w
end

# ----------------------------------------------------------------------------------
# Design write-back (set_design! / constructor)
# ----------------------------------------------------------------------------------

"""
    set_design!(sm; Ku=..., Aex=..., D=...)

Write new design values into the underlying `sim` and the `StateMatch`'s design
declaration.  Scalars/`Fill`s keep their O(1) storage, arrays/functions are
materialised through `make_param` (the same path `add_exch`/`add_anis`/`add_dmi`
use).

Stencil-cache invalidation on write-back: the exchange/DMI pair tables in the
stencil layer key only on `mesh.layout_version`, so a parameter write must drop
them explicitly — `set_design!` clears the pair tables (and the interfacial
`Dcls`) of the interaction it rewrites, mirroring what `set_Ms` does via
`_drop_ms_caches!` but without discarding the (Ms-only) class map.  Uniaxial
anisotropy is a pointwise kernel with no caches.
"""
function set_design!(sm::StateMatch; kwargs...)
    for k in keys(kwargs)
        k in _STATEMATCH_KEYS || throw(ArgumentError(
            "unsupported design key :$k (supported: $(_STATEMATCH_KEYS))"))
    end
    sm.design = merge(sm.design, kwargs)
    _apply_design!(sm)
    return sm
end

function _apply_design!(sm::StateMatch)
    sim = sm.sim
    T = eltype(sim.spin)
    if :Ku in keys(sm.design)
        sm.anis = sm.anis === nothing ? _find_interaction(Anisotropy, sim, :Ku) : sm.anis
        _write_Ku!(sm.anis, sm.design.Ku, sim, T)
    end
    if :Aex in keys(sm.design)
        sm.exch = sm.exch === nothing ? _find_interaction(Exchange, sim, :Aex) : sm.exch
        _write_Aex!(sm.exch, sm.design.Aex, sim, T)
    end
    if :D in keys(sm.design)
        sm.dmi = sm.dmi === nothing ? _find_interaction(DMI, sim, :D) : sm.dmi
        _write_D!(sm.dmi, sm.design.D, sim, T)
    end
    return sm
end

# Ku: pointwise kernel reads anis.Ku every call — no caches to drop.
_write_Ku!(anis::Anisotropy, v, sim::AbstractSim, T) =
    (anis.Ku = make_param(T, v, sim.mesh, sim.n_total); anis)

# Aex: the scalar-per-cell convention feeds Ax = Ay = Az; the pair tables built
# by the stencil layer key on mesh.layout_version only, so drop them explicitly.
function _write_Aex!(exch::Exchange, v, sim::AbstractSim, T)
    a = make_param(T, v, sim.mesh, sim.n_total)
    exch.Ax = exch.Ay = exch.Az = a
    exch.pair_Ax = exch.pair_Ay = exch.pair_Az = nothing
    exch.stencil_layout = -1
    return exch
end

# D: one scalar D per cell feeding Dx = Dy = Dz (the `add_dmi(sim, D::Number)`
# convention; the grad_dmi kernels assume it).  Same pair-table drop as Aex,
# plus the interfacial per-class D guard.
function _write_D!(dmi::DMI, v, sim::AbstractSim, T)
    d = make_param(T, v, sim.mesh, sim.n_total)
    dmi.Dx = dmi.Dy = dmi.Dz = d
    dmi.pair_Dx = dmi.pair_Dy = dmi.pair_Dz = nothing
    empty!(dmi.Dcls)
    dmi.stencil_layout = -1
    return dmi
end

# ----------------------------------------------------------------------------------
# Forward
# ----------------------------------------------------------------------------------

"""
    run_forward!(sm) -> J

Forward pass: write the current design into `sim`, reset to the fixed initial
state, relax to a steady state, and freeze `m*`, the total field `H0` and the
adjoint workspace for [`gradient!`](@ref).  Returns the loss
`J = Σ w_i |m*_i − target_i|²` (vacuum cells excluded).

The relaxation residual is recorded as `sm.max_transverse` (max per-site
|m*×H0|, the transverse field that should vanish at equilibrium); a warning is
printed when it is too large for trustworthy gradients (> 1e-6 of the field
scale) — tighten `stopping_dmdt` in that case.
"""
function run_forward!(sm::StateMatch)
    sim = sm.sim
    # Deterministic forwards: the adaptive integrator carries step-size and
    # FSAL state across relax calls — a stale k7 poisons the first step after
    # init_m0 resets y_current, and a run-away step_next can stop a fresh
    # relax instantly.  Reset both so every forward is the same fresh
    # trajectory from the fixed initial state (a StateMatch loss J(θ) is only
    # well-defined if the forward is history-independent).
    if hasproperty(sim.driver, :integrator)
        integ = sim.driver.integrator
        integ.step_next = 0.0
        integ.succeed = false
        integ.t = 0.0
    end
    _apply_design!(sm)                       # idempotent; keeps FD perturbations honest
    init_m0(sim, sm.init_m)
    relax(sim; max_steps = sm.max_steps, stopping_dmdt = sm.stopping_dmdt)

    sm.m_star = Float64.(Array(sim.spin))
    effective_field(sim, sim.spin, 0.0)      # total field at m* (all interactions)
    sm.H0 = Float64.(Array(sim.field))
    sm.ws = AdjointWorkspace(sim)            # Ms tiles must match the current design

    # convergence diagnostic: the transverse field vanishes at a steady state
    N = sim.n_total
    hmax = 0.0
    mtrans = 0.0
    @inbounds for i in 1:N
        j = 3i - 2
        hx, hy, hz = sm.H0[j], sm.H0[j+1], sm.H0[j+2]
        hn = abs(hx) + abs(hy) + abs(hz)
        hn > hmax && (hmax = hn)
        mx, my, mz = sm.m_star[j], sm.m_star[j+1], sm.m_star[j+2]
        tx = my * hz - mz * hy
        ty = mz * hx - mx * hz
        tz = mx * hy - my * hx
        tn = sqrt(tx^2 + ty^2 + tz^2)
        tn > mtrans && (mtrans = tn)
    end
    sm.max_transverse = mtrans
    if mtrans > 1e-6 * hmax
        @warn "StateMatch: the relaxed state is far from a steady state " *
              "(max |m×H| = $(@sprintf("%.3e", mtrans)) vs field scale " *
              "$(@sprintf("%.3e", hmax))); the gradient! result will be " *
              "inaccurate. Tighten stopping_dmdt."
    end

    sm.J = _loss_value(sm)
    return sm.J
end

function _loss_value(sm::StateMatch)
    N = sm.sim.n_total
    w = _loss_weights(sm)
    m, t = sm.m_star, sm.target_m
    J = 0.0
    @inbounds for i in 1:N
        j = 3i - 2
        J += w[i] * (abs2(m[j] - t[j]) + abs2(m[j+1] - t[j+1]) + abs2(m[j+2] - t[j+2]))
    end
    return J
end

# ----------------------------------------------------------------------------------
# Backward
# ----------------------------------------------------------------------------------

"""
    gradient!(sm) -> G

Backward pass (requires [`run_forward!`](@ref) first): one steady-state adjoint
solve plus the ∂f/∂θ accumulation.  Returns a `NamedTuple` shaped like the
declared design — a scalar for scalar/`Fill` design values, the per-cell map
(length `n_total`) for array/function values.  The result is also cached in
`sm.grad`.

Conventions (validated against central finite differences in
`test/adjoint/test_statematch_fd.jl`): with `J = Σ w‖m* − target‖²` the loss
gradient fed to the adjoint solve is `2·w·P_t(m* − target)`; the design
gradient is `∂J/∂θ = w·(∂H/∂θ)` with the reduced weight
`w = γ′(λ×m* − α m*×(λ×m*))`, i.e. the −λᵀ∂f/∂θ sign and the γ′ factor are
folded into `w`.
"""
function gradient!(sm::StateMatch)
    sm.m_star !== nothing || error("run_forward! must be called before gradient!")
    sim = sm.sim
    N = sim.n_total
    n3 = 3N
    m_star, H0, ws = sm.m_star, sm.H0, sm.ws
    back = get_backend(sim.spin)
    alpha, gamma = _adjoint_alpha_gamma(sim)

    # 1. loss gradient: kernel gives P_t(m* − target); J = Σ w‖·‖² → scale by 2·w
    grad_g = zeros(n3)
    MicroMagnetic.loss_grad_kernel!(back, groupsize[])(grad_g, m_star, sm.target_m;
                                                       ndrange = N)
    w3 = repeat(2 .* _loss_weights(sm); inner = 3)
    grad_g .*= w3

    # 2. steady-state adjoint: (Df)ᵀλ = P_t ∇g at the frozen m*
    vjp! = (out, lam) -> llg_rhs_vjp!(out, sim, m_star, H0, lam, ws;
                                      alpha = alpha, gamma = gamma)
    lam = solve_steady_adjoint(vjp!, grad_g, m_star; backend = sm.adjoint_backend,
                               tol = sm.adjoint_tol, maxiter = sm.adjoint_maxiter,
                               verbose = sm.verbose)

    # 3. reduced θ-weight (sign & γ′ folded in)
    w = zeros(n3)
    gamma_prime = gamma / (1 + alpha^2)
    MicroMagnetic.llg_param_weight_kernel!(back, groupsize[])(w, m_star, lam, alpha,
                                                              gamma_prime; ndrange = N)

    # 4. per-key ∂J/∂θ accumulation
    Gs = map(keys(sm.design)) do key
        g = _grad_key!(zeros(N), sm, key, w, back)
        v = getfield(sm.design, key)
        v isa Number || v isa Fill ? sum(g) : g
    end
    sm.grad = NamedTuple{keys(sm.design)}(Tuple(Gs))
    return sm.grad
end

function _grad_key!(G::Vector{Float64}, sm::StateMatch, key::Symbol, w::Vector{Float64},
                    back)
    sim = sm.sim
    N = sim.n_total
    if key === :Ku
        anis = sm.anis
        MicroMagnetic.grad_ku_kernel!(back, groupsize[])(G, sm.m_star, w, anis.axis_x,
                                                         anis.axis_y, anis.axis_z,
                                                         sim.mu0_Ms; ndrange = N)
    elseif key === :Aex
        exch = sm.exch
        mesh = sim.mesh
        MicroMagnetic.grad_aex_kernel!(back, _stencil_wg(mesh))(G, sm.m_star, w,
            sim.mu0_Ms, sim.inv_ms, exch.Ax, exch.Ay, exch.Az,
            Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
            mesh.nx, mesh.ny, mesh.nz, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
            ndrange = (mesh.nx, mesh.ny, mesh.nz))
    elseif key === :D
        dmi = sm.dmi
        mesh = sim.mesh
        if dmi.type === :bulk
            MicroMagnetic.grad_dmi_kernel!(back, _stencil_wg(mesh))(G, sm.m_star, w,
                sim.mu0_Ms, sim.inv_ms, dmi.Dx, dmi.Dy, dmi.Dz,
                Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                mesh.nx, mesh.ny, mesh.nz, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                ndrange = (mesh.nx, mesh.ny, mesh.nz))
        else
            MicroMagnetic.grad_dmi_interfacial_kernel!(back, _stencil_wg(mesh))(G,
                sm.m_star, w, sim.mu0_Ms, sim.inv_ms, dmi.Dx,
                Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                mesh.nx, mesh.ny, mesh.nz, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                ndrange = (mesh.nx, mesh.ny, mesh.nz))
        end
    else
        throw(ArgumentError("unsupported design key :$key"))
    end
    return G
end

"""
    value_and_gradient(sm) -> (J, G)

Convenience pair: `run_forward!(sm)` followed by `gradient!(sm)`.
"""
function value_and_gradient(sm::StateMatch)
    J = run_forward!(sm)
    G = gradient!(sm)
    return (J, G)
end
