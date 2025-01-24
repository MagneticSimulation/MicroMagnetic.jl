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

B = build_matrix(sim, gamma=2.211e5)

function compute_eigen_values(matrix, sim; nev=10)
    sparse_matrix = sparse(matrix)
    vals, vecs = Arpack.eigs(sparse_matrix, nev=2*nev, which=:SM, tol=1e-6, maxiter=10*sim.n_total)

    indices = findall(x -> imag(x) > 0, vals)
    freqs = vals[indices]
    eigenvectors = vecs[:, indices]

    mesh = sim.mesh
    N = sim.n_total
    evecs = reshape(eigenvectors, (2, N, nev))
    new_evecs = zeros(eltype(evecs), 3, N, nev)
    new_evecs[1:2, :, :] .= evecs

    m0 = reshape(sim.spin, 3, N)
    m = zeros(Complex{Float64}, 3, N, nev)
    for i = 1:N
        R = MicroMagnetic.rotation_matrix(m0[1, i], m0[2, i], m0[3, i])
        for j = 1:nev
            m[:, i, j] = R * new_evecs[:, i, j]
        end
    end
    
    m = reshape(m, (3, mesh.nx, mesh.ny, mesh.nz, nev))

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

freqs, m = compute_eigen_values(B, sim, nev=10)
plot(freqs, m)
