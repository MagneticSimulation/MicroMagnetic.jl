# Acceptance tests for the steady-state adjoint solver (src/adjoint/steady.jl).
#
# Toy operators mimic the LLG linearization structure:
#
#     A = P_t [ −c (I + βL) + (C L + L C)/2 ] P_t ,   C = blockdiag([m̄ᵢ]×)
#
# with L a symmetric exchange-like Laplacian (1D ring or 2D 5-point stencil),
# c > 0 the damping strength (symmetric negative-definite part), the
# symmetrized cross-product term the antisymmetric (precession) part, and P_t
# the per-site tangential projector.  The outer P_t gives A the exact
# per-site null space of the true LLG adjoint at equilibrium:
# A m̄ = Aᵀ m̄ = 0.  The operator is dense here (N ≤ 64), but the solver only
# ever sees it through the matrix-free VJP closure.
using MicroMagnetic
using Test
using LinearAlgebra
using Random

_cross_matrix(v) = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]

# symmetric ring Laplacian (exchange-like stencil) of strength s
function _toy_laplacian_1d(N::Int, s::Float64)
    L = zeros(N, N)
    for i in 1:N
        L[i, i] += 2s
        L[i, mod1(i - 1, N)] -= s
        L[i, mod1(i + 1, N)] -= s
    end
    return L
end

# 5-point stencil on a periodic nx×nx grid
function _toy_laplacian_2d(nx::Int, s::Float64)
    L = zeros(nx * nx, nx * nx)
    idx(i, j) = (i - 1) * nx + j
    for i in 1:nx, j in 1:nx
        L[idx(i, j), idx(i, j)] += 4s
        L[idx(i, j), idx(mod1(i - 1, nx), j)] -= s
        L[idx(i, j), idx(mod1(i + 1, nx), j)] -= s
        L[idx(i, j), idx(i, mod1(j - 1, nx))] -= s
        L[idx(i, j), idx(i, mod1(j + 1, nx))] -= s
    end
    return L
end

function _tangential_projector(m::AbstractVector{Float64})
    P = zeros(length(m), length(m))
    for i in 1:3:length(m)
        for a in 1:3, b in 1:3
            P[i + a - 1, i + b - 1] = (a == b ? 1.0 : 0.0) - m[i + a - 1] * m[i + b - 1]
        end
    end
    return P
end

function _random_unit_field!(rng, m::AbstractVector{Float64})
    randn!(rng, m)
    for i in 1:3:length(m)
        m[i:i+2] ./= norm(view(m, i:i+2))
    end
    return m
end

"""
    _toy_adjoint_matrix((:ring, N) | (:grid, nx); seed, c, beta, ls) -> (A, m)

Dense toy adjoint operator with LLG structure: symmetric negative-definite
damping `−c(I + βL)`, antisymmetric cross-product part `(CL + LC)/2`, per-site
null space along `m̄`.  `(CL + LC)/2` is the antisymmetrization of the
precession-adjoint `C·L` and is exactly antisymmetric since
`([m̄]×L + L[m̄]×)ᵀ = −([m̄]×L + L[m̄]×)`.
"""
function _toy_adjoint_matrix(lattice::Tuple{Symbol,Int}; seed::Int = 20260905,
                             c::Float64 = 2.0, beta::Float64 = 1.0, ls::Float64 = 0.5)
    kind, N = lattice
    N >= 2 || throw(ArgumentError("need at least 2 sites"))
    L = kind === :ring ? _toy_laplacian_1d(N, ls) : _toy_laplacian_2d(isqrt(N), ls)
    n3 = 3N
    rng = MersenneTwister(seed)
    m = _random_unit_field!(rng, Vector{Float64}(undef, n3))
    L3 = kron(Matrix{Float64}(I(3)), L)
    C = zeros(n3, n3)
    for i in 1:3:n3
        C[i:i+2, i:i+2] .= _cross_matrix(view(m, i:i+2))
    end
    A = _tangential_projector(m) *
        (-c * (Matrix{Float64}(I, n3, n3) + beta * L3) + (C * L3 + L3 * C) / 2) *
        _tangential_projector(m)
    return A, m
end

_toy_vjp(A) = (out, x) -> mul!(out, A, x)

# largest |λᵢ · mᵢ| over sites (tangency defect)
_max_site_dot(λ, m) =
    maximum(abs(dot(view(λ, i:i+2), view(m, i:i+2))) for i in 1:3:length(m))

@testset "adjoint/steady" begin

    @testset "tangential projection" begin
        rng = MersenneTwister(11)
        N = 10
        m = _random_unit_field!(rng, Vector{Float64}(undef, 3N))
        x = randn(rng, 3N)
        px = project_tangent(x, m)
        @test px ≈ _tangential_projector(m) * x
        @test _max_site_dot(px, m) < 1e-14
        @test norm(project_tangent(3.0 * m, m), Inf) < 1e-14  # longitudinal annihilated
        y = copy(x)
        project_tangent!(y, y, m)                        # out may alias x
        @test y == px
    end

    @testset "toy operator sanity" begin
        A, m = _toy_adjoint_matrix((:ring, 16))
        @test norm(A * m) < 1e-12                        # per-site null space, both sides
        @test norm(A' * m) < 1e-12
        @test maximum(real.(eigvals(A))) < 1e-12         # no growth (left half-plane)
        rng = MersenneTwister(7)
        x = _tangential_projector(m) * randn(rng, length(m))
        @test dot(x, A * x) < 0.0                        # dissipative quadratic form
    end

    @testset "residual, tangency, backend agreement ($kind N=$N)" for (kind, N) in
        [(:ring, 2), (:ring, 5), (:ring, 16), (:ring, 64), (:grid, 64)]
        A, m = _toy_adjoint_matrix((kind, N); seed = 20260905 + N)
        vjp! = _toy_vjp(A)
        rng = MersenneTwister(99 + N)
        g = randn(rng, 3N)
        tol = 1e-8

        λpt = solve_steady_adjoint(vjp!, g, m; backend = :pseudotime, tol)
        λkr = solve_steady_adjoint(vjp!, g, m; backend = :krylov, tol)

        b = project_tangent(g, m)
        for λ in (λpt, λkr)
            r = similar(g)
            vjp!(r, λ)
            r .-= b
            @test norm(r) <= tol * norm(g)                               # (a) residual
            @test _max_site_dot(λ, m) <= 1e-12 * max(1.0, norm(λ, Inf))  # (b) tangency
        end
        @test norm(λpt - λkr) / max(norm(λpt), norm(λkr)) < 1e-6         # (c) agreement
    end

    @testset "inversion recovers λ_true (N=$N)" for N in (16, 64)
        A, m = _toy_adjoint_matrix((:ring, N); seed = 4242 + N)
        vjp! = _toy_vjp(A)
        rng = MersenneTwister(555 + N)
        λtrue = project_tangent(randn(rng, 3N), m)       # tangential ground truth
        g = similar(λtrue)
        vjp!(g, λtrue)                                   # ∇g = (Df)ᵀ λ_true
        tol = 1e-10
        λpt = solve_steady_adjoint(vjp!, g, m; backend = :pseudotime, tol)
        λkr = solve_steady_adjoint(vjp!, g, m; backend = :krylov, tol)
        @test norm(λpt - λtrue) / norm(λtrue) < 1e-6
        @test norm(λkr - λtrue) / norm(λtrue) < 1e-6
    end

    @testset "edge cases" begin
        A, m = _toy_adjoint_matrix((:ring, 8); seed = 31)
        vjp! = _toy_vjp(A)
        n3 = length(m)
        g = randn(MersenneTwister(5), n3)

        # zero gradient / purely longitudinal gradient → zero adjoint
        @test solve_steady_adjoint(vjp!, zeros(n3), m) == zeros(n3)
        @test solve_steady_adjoint(vjp!, zeros(n3), m; backend = :krylov) == zeros(n3)
        @test solve_steady_adjoint(vjp!, 2.0 * m, m) == zeros(n3)
        @test solve_steady_adjoint(vjp!, 2.0 * m, m; backend = :krylov) == zeros(n3)

        # malformed input
        @test_throws DimensionMismatch solve_steady_adjoint(vjp!, zeros(n3 + 1), m)
        @test_throws ArgumentError solve_steady_adjoint(vjp!, g, 2.0 .* m)  # non-unit m
        # zero (vacuum) sites are legal — with a gradient consistent
        # with them (vanishing vacuum rows) the solve succeeds and the vacuum
        # adjoint components stay zero
        mvac = copy(m)
        mvac[1:3] .= 0.0
        gvac = 2.0 .* mvac                       # purely longitudinal, 0 at vacuum
        @test solve_steady_adjoint(vjp!, gvac, mvac) == zeros(n3)
        @test solve_steady_adjoint(vjp!, gvac, mvac; backend = :krylov) == zeros(n3)
        @test_throws ArgumentError solve_steady_adjoint(vjp!, g, m; backend = :nope)
        @test_throws ErrorException solve_steady_adjoint(vjp!, g, m;
                                                         backend = :pseudotime, maxiter = 1)
    end
end
