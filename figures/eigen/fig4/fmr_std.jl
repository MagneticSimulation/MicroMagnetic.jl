using MicroMagnetic

#using LinearAlgebra
using SparseArrays
using Arpack
using CairoMakie
using Printf
#@using_gpu()
function setup(;m0=(0,0,1))
    mesh = FDMesh(nx=24, ny=24, nz=2, dx=5e-9, dy=5e-9, dz=5e-9)
    v = (1, 0.7, 0)
    H = 80e3 .* v ./ sqrt(sum(v.^2))
    sim = create_sim(mesh; A=1.3e-11, H=H, m0=m0, demag=true, Ms=8e5)
    return sim
end

# we used Float64 to compute the static state
MicroMagnetic.set_precision(Float64)
sim = setup()
relax(sim, stopping_dmdt=0.0001)

m0 = Array(sim.spin)

MicroMagnetic.set_precision(AbstractFloat)
sim = setup(m0=m0)

# Build phase: matrix-free operator (O(N) storage — no 2N×2N dense Jacobian).
# The returned LLGJacOperator supports size/eltype, mul!(3/5-arg), issymmetric
# and copy(), so it can be passed directly to Arpack.eigs for Krylov-targeted
# eigensolves (which=:LR / :LM / :SR).  When targeting interior eigenvalues via
# shift-invert (e.g. which=:SM → smallest-magnitude lowest-freq FMR modes)
# Arpack requires a factorisation of (B−σI), so we materialise a sparse matrix
# once via the built-in SparseArrays.sparse(op) (= 2N mul! calls).  For truly
# large meshes (2N ≫ 10⁴) and :LR / :LM targeting, call
#     Arpack.eigs(op; nev=..., which=:LR, explicittransform=:none)
# directly on the operator and skip the sparse() step.
op = build_matrix(sim, gamma=2.211e5, matrixfree=true)

function compute_eigen_values(op, sim; nev=10, which=:SM)
    if which == :SM
        # Smallest-magnitude (lowest-freq FMR) modes live in the interior of
        # the spectrum; Arpack needs shift-invert → factorise via sparse.
        A = SparseArrays.sparse(op)
        vals, vecs = Arpack.eigs(A; nev=2*nev, which=:SM, tol=1e-6,
                                 maxiter=10*sim.n_total, explicittransform=:none)
    else
        # Pure matrix-free Krylov path (e.g. :LR → least-damped modes).
        vals, vecs = Arpack.eigs(op; nev=2*nev, which=which, tol=1e-6,
                                 maxiter=10*sim.n_total, explicittransform=:none)
    end

    indices = findall(x -> imag(x) > 0, vals)
    freqs = vals[indices]
    eigenvectors = vecs[:, indices]
    # Ensure we got at least nev valid positive-imaginary eigenvalues.
    if length(freqs) < nev
        @warn("Only $(length(freqs))/$(nev) positive-imaginary modes returned " *
              "by Arpack; increase nev or maxiter.")
        keep = length(freqs)
    else
        keep = nev
    end
    freqs = freqs[1:keep]
    eigenvectors = eigenvectors[:, 1:keep]

    mesh = sim.mesh
    N = sim.n_total
    evecs = reshape(eigenvectors, (2, N, keep))
    new_evecs = zeros(eltype(evecs), 3, N, keep)
    new_evecs[1:2, :, :] .= evecs

    m0 = reshape(sim.spin, 3, N)
    m = zeros(Complex{Float64}, 3, N, keep)
    for i = 1:N
        R = MicroMagnetic.rotation_matrix(m0[1, i], m0[2, i], m0[3, i])
        for j = 1:keep
            m[:, i, j] = R * new_evecs[:, i, j]
        end
    end
    
    m = reshape(m, (3, mesh.nx, mesh.ny, mesh.nz, keep))

    return freqs, m
end

function plot(freqs, m)

    frequencies_ghz = imag(freqs) ./ (2pi * 1e9)                                

    fig = Figure(size = (600, 290), fontsize=14)

    ncols = 5  
    nrows = ceil(Int, length(freqs) / ncols)  

    for i in 1:length(freqs)
        vv = abs.(m[2, :, :, 1, i]) 
        fr = frequencies_ghz[i]
        freq = @sprintf("%.3f GHz", fr)
        ax = Axis(fig[div(i-1, ncols)+1, (i-1) % ncols + 1], title = " $(freq)", aspect = 1)
        hidedecorations!(ax)
        heatmap!(ax, vv, colormap = :coolwarm, interpolate=true)
    end

    save("fig4.pdf", fig)
end

freqs, m = compute_eigen_values(op, sim, nev=10, which=:SM)
plot(freqs, m)
