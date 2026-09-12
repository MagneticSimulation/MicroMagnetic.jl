# ---------------------------------------------------------------------------
# Steady-state adjoint solver.
#
# At a frozen steady state m* (f(m*) = 0, per-site unit vectors) the stage
# gradient of an inverse problem is the solution of one linear system
#
#     (Df)ᵀ λ = P_t ∇g ,     P_t = I − m* m*ᵀ (per-site tangential projector)
#
# LLG preserves |m|, so the stacked per-site state m̄ is a null vector of both
# Df and (Df)ᵀ at equilibrium: the system is singular along the per-site
# longitudinal directions.  ∇g is therefore projected tangentially first and
# λ is kept tangential (initial guess 0); the tangential solution is unique
# whenever the symmetric part of (Df)ᵀ is negative definite on the tangent
# space (the α > 0 damping analogue).
#
# The VJP is supplied by the caller: vjp!(out, x) computes out = (Df)ᵀx at
# the frozen m* (typically linearize.jl's apply_KT!).  This solver is
# agnostic to how the VJP is evaluated and never materializes a matrix.
#
# Backends:
#   * :pseudotime — explicit adaptive pseudo-time marching of the flow
#     λ̇ = (Df)ᵀλ − P_t∇g (Richardson iteration), stopped when the residual
#     ‖(Df)ᵀλ − P_t∇g‖ ≤ tol·‖P_t∇g‖.  Sign note: spec((Df)ᵀ) = spec(Df) lies
#     in the left half-plane for a stable forward state, so *this* is the
#     convergent direction of pseudo-time; the flow with the opposite sign
#     (λ̇ = −(Df)ᵀλ + P_t∇g) reaches the same steady state only when
#     integrated in reversed pseudo-time.
#   * :krylov — GMRES (Krylov.jl) acting on the matrix-free operator
#     P_t (Df)ᵀ P_t, with the same relative residual tolerance.
#
# Float64/CPU only at this stage; the KA/GPU generalization follows with the
# kernel layer (T1).
# ---------------------------------------------------------------------------

using LinearAlgebra
using Krylov

"""
    project_tangent!(out, x, m_star)
    project_tangent(x, m_star)

Per-site tangential projection `P_t x` with the projector
`P_t = I − m* m*ᵀ`: for each stacked site `i`,

    out[3i-2:3i] = x[3i-2:3i] − (x[3i-2:3i] · m_star[3i-2:3i]) m_star[3i-2:3i]

`out` may alias `x`.  Vacuum sites (`m = 0`) pass through unchanged
(`P_t = I` there).
"""
function project_tangent!(out::AbstractVector{Float64}, x::AbstractVector{Float64},
                          m_star::AbstractVector{Float64})
    @inbounds for i in 1:3:length(m_star)
        mx = m_star[i]
        my = m_star[i + 1]
        mz = m_star[i + 2]
        dx = x[i] * mx + x[i + 1] * my + x[i + 2] * mz
        out[i] = x[i] - dx * mx
        out[i + 1] = x[i + 1] - dx * my
        out[i + 2] = x[i + 2] - dx * mz
    end
    return out
end

project_tangent(x::AbstractVector{Float64}, m_star::AbstractVector{Float64}) =
    project_tangent!(Vector{Float64}(undef, length(x)), x, m_star)

function _validate_m_star(m_star::Vector{Float64})
    for (site, i) in enumerate(1:3:length(m_star))
        nrm = sqrt(abs2(m_star[i]) + abs2(m_star[i + 1]) + abs2(m_star[i + 2]))
        # vacuum sites (|m| = 0, e.g. Ms = 0 cells in a StateMatch problem) are
        # legitimate: their (Df)ᵀ rows vanish and P_t passes them through.
        if nrm > 1e-6 && abs(nrm - 1.0) > 1e-6
            throw(ArgumentError("m_star must hold unit vectors per site (vacuum " *
                                "sites may be zero); " *
                                "site $site has |m| = $nrm"))
        end
    end
    return nothing
end

# Matrix-free operator for Krylov.jl: mul! is the only supported operation.
struct _TangentAdjointOp{F} <: AbstractMatrix{Float64}
    n::Int
    apply!::F
end

Base.size(op::_TangentAdjointOp) = (op.n, op.n)
Base.getindex(op::_TangentAdjointOp, i::Int, j::Int) =
    error("_TangentAdjointOp is matrix-free; only mul! is supported")

function LinearAlgebra.mul!(y::AbstractVector{Float64}, op::_TangentAdjointOp,
                            x::AbstractVector{Float64})
    return op.apply!(y, x)
end

"""
    solve_steady_adjoint(vjp!, grad_g, m_star; backend=:pseudotime, tol=1e-8,
                         maxiter=10_000, dt0=0.0, memory=0, verbose=false)
                         -> Vector{Float64}

Solve the steady-state adjoint equation

    (Df)ᵀ λ = P_t ∇g ,     P_t = I − m* m*ᵀ (per site)

for the adjoint field λ at the frozen steady state `m*` (`f(m*) = 0`,
stacked 3N unit vectors; vacuum sites with `m = 0` are allowed — their
`(Df)ᵀ` rows and `P_t`-projected gradient must vanish there, which holds
whenever the loss excludes vacuum cells).

`vjp!(out, x)` must compute `out = (Df)ᵀ x` at the frozen `m*` (a caller
closure; typically linearize.jl's apply_KT!).  The solver never materializes
`(Df)ᵀ`.

Singularity handling: LLG preserves |m|, hence `m*` spans a per-site null
space of `(Df)ᵀ` at equilibrium.  `grad_g` is projected onto the tangent
space and λ is kept tangential; when the symmetric part of `(Df)ᵀ` is
negative definite on the tangent space (α > 0), the tangential solution is
unique and both backends converge to it.

Backends:
- `:pseudotime` (default): explicit adaptive pseudo-time marching of
  `λ̇ = (Df)ᵀλ − P_t∇g` with the residual projected onto the tangent space each
  step (with a non-homogeneous field term such as Zeeman the longitudinal
  directions are not in the kernel of `(Df)ᵀ`, and the unprojected residual
  norm would never contract); stops when the tangential residual
  `‖P_t((Df)ᵀλ − P_t∇g)‖ ≤ tol·‖P_t∇g‖`.
- `:krylov`: GMRES (Krylov.jl) on the matrix-free operator `P_t (Df)ᵀ P_t`,
  same relative residual tolerance.

Both backends return the same tangential solution up to O(tol) relative
error.  The longitudinal components of `grad_g` never influence the result.

# Extended kwargs
- `dt0`: initial pseudo-time step; `0` (default) estimates it from one VJP
  probe and an adaptive reject/halve loop stabilizes the march.
- `memory`: GMRES restart length; `0` (default) uses `min(3N, 100)`.
- `verbose`: log backend, iteration count and final relative residual.

Throws `ArgumentError` for malformed input and `ErrorException` when the
selected backend fails to reach `tol` within `maxiter` iterations.
"""
function solve_steady_adjoint(vjp!::Function, grad_g::Vector{Float64},
                              m_star::Vector{Float64};
                              backend::Symbol = :pseudotime,
                              tol::Float64 = 1e-8,
                              maxiter::Int = 10_000,
                              dt0::Float64 = 0.0,
                              memory::Int = 0,
                              verbose::Bool = false)
    length(grad_g) == length(m_star) ||
        throw(DimensionMismatch("grad_g and m_star must have equal length " *
                                "($(length(grad_g)) vs $(length(m_star)))"))
    n3 = length(grad_g)
    n3 % 3 == 0 || throw(ArgumentError("grad_g must be a stacked per-site vector field " *
                                       "(length divisible by 3), got $n3"))
    backend in (:pseudotime, :krylov) ||
        throw(ArgumentError("backend must be :pseudotime or :krylov (got $backend)"))
    tol > 0 || throw(ArgumentError("tol must be positive (got $tol)"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive (got $maxiter)"))
    _validate_m_star(m_star)

    b = project_tangent(grad_g, m_star)                 # P_t ∇g
    bnorm = norm(b)
    # A numerically longitudinal ∇g leaves only projection roundoff in b; the
    # adjoint field is then 0 up to that roundoff, and no backend can reach a
    # relative tolerance against such a right-hand side.
    if bnorm <= 128 * eps(Float64) * max(norm(grad_g), 1.0)
        verbose && @info "solve_steady_adjoint: P_t ∇g ≈ 0 (numerically longitudinal); λ = 0."
        return zeros(n3)
    end

    tstart = time()
    if backend === :pseudotime
        λ, iters, relres = _pseudotime_adjoint(vjp!, b, m_star, tol; maxiter, dt0)
    else
        λ, iters, relres = _krylov_adjoint(vjp!, b, m_star, tol; maxiter, memory)
    end
    project_tangent!(λ, λ, m_star)                      # strip roundoff drift
    verbose && @info("solve_steady_adjoint finished", backend = backend,
                     iterations = iters, relative_residual = relres,
                     seconds = round(time() - tstart; digits = 3))
    return λ
end

# λ̇ = (Df)ᵀλ − b integrated with an explicit Euler/Richardson scheme:
#     λₙ₊₁ = λₙ + dt·rₙ ,   rₙ = P_t((Df)ᵀλₙ − b),
# with a residual-controlled step: a step is accepted when the residual does
# not grow, otherwise dt is halved and the step is retried.  For a dissipative
# operator (negative-definite symmetric part on the tangent space) a
# contracting dt always exists, since
#     ‖rₙ₊₁‖² = ‖rₙ‖² + 2dt·rₙᵀA rₙ + dt²‖A rₙ‖²  and  rₙᵀA rₙ = rₙᵀS rₙ < 0.
#
# The per-step tangential projection is load-bearing, not cosmetic: with a
# non-homogeneous field term (e.g. Zeeman) the radial direction is not in the
# kernel of Df (f(c·m) ≠ 0 for c ≠ 1), so (Df)ᵀ leaks tangential iterates into
# the longitudinal directions and the unprojected residual norm never
# contracts.  Projecting keeps the march on the tangent space, where the
# damped-LLG contraction argument holds; the resulting fixed point
# P_t((Df)ᵀλ − b) = 0 is exactly the equation the :krylov backend solves.
function _pseudotime_adjoint(vjp!::Function, b::Vector{Float64},
                             m_star::Vector{Float64}, tol::Float64;
                             maxiter::Int, dt0::Float64)
    n = length(b)
    λ = zeros(n)
    r = Vector{Float64}(undef, n)          # current residual
    λtrial = Vector{Float64}(undef, n)
    rtrial = Vector{Float64}(undef, n)

    bnorm = norm(b)
    stop = tol * bnorm
    @. r = -b                              # residual of the guess λ = 0
    rnorm = bnorm

    dt = if dt0 > 0.0
        dt0
    else
        @. rtrial = b / bnorm              # probe the operator scale
        vjp!(λtrial, rtrial)
        est = norm(λtrial)
        isfinite(est) || throw(ErrorException("vjp! produced non-finite values"))
        0.5 / max(est, eps(Float64))
    end
    dt_start = dt

    iters = 0
    while iters < maxiter
        rnorm <= stop && break
        @. λtrial = λ + dt * r
        vjp!(rtrial, λtrial)
        rtrial .-= b
        project_tangent!(rtrial, rtrial, m_star)
        iters += 1
        ntrial = norm(rtrial)
        isfinite(ntrial) ||
            throw(ErrorException("vjp! produced non-finite values during the pseudo-time march"))
        if ntrial <= rnorm
            λ, λtrial = λtrial, λ          # accept the step
            r, rtrial = rtrial, r
            rnorm = ntrial
            dt = min(1.2 * dt, 1e8 * dt_start)
        else
            dt *= 0.5
            dt <= 1e-6 * dt_start && throw(ErrorException(
                "solve_steady_adjoint (:pseudotime): the pseudo-time step collapsed " *
                "(dt = $dt). Is the symmetric part of (Df)ᵀ negative definite on the " *
                "tangent space (the α > 0 analogue)?"))
        end
    end
    rnorm <= stop || throw(ErrorException(
        "solve_steady_adjoint (:pseudotime) did not converge in $maxiter iterations " *
        "(relative residual $(rnorm / bnorm), tol $tol)"))
    return λ, iters, rnorm / bnorm
end

# GMRES on the tangentially projected operator P_t (Df)ᵀ P_t.  The operator
# maps everything into the tangent space, so the Krylov space — and hence the
# iterate — stays tangential for a tangential right-hand side.
function _krylov_adjoint(vjp!::Function, b::Vector{Float64}, m_star::Vector{Float64},
                         tol::Float64; maxiter::Int, memory::Int)
    n = length(b)
    work = Vector{Float64}(undef, n)
    apply! = function (y, x)
        project_tangent!(work, x, m_star)  # work = P_t x
        vjp!(y, work)                      # y = (Df)ᵀ P_t x
        project_tangent!(y, y, m_star)     # enforce exact tangency
        return y
    end
    op = _TangentAdjointOp(n, apply!)
    mem = memory > 0 ? min(memory, n) : min(n, 100)
    λ, stats = Krylov.gmres(op, b; atol = 0.0, rtol = tol, memory = mem, itmax = maxiter)
    stats.solved || throw(ErrorException("solve_steady_adjoint (:krylov) failed: " *
                                         string(stats.status)))
    # Independent residual of the *unprojected* equation for reporting/tests.
    res = Vector{Float64}(undef, n)
    vjp!(res, λ)
    res .-= b
    return λ, Int(stats.niter), norm(res) / norm(b)
end

export solve_steady_adjoint, project_tangent, project_tangent!
