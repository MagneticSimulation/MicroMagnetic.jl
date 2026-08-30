```@meta
ShareDefaultModule = true
```

# Eigenmode Analysis

MicroMagnetic.jl provides tools for computing the eigenmodes of a magnetic
system around a static equilibrium state. The linearised Landau–Lifshitz–Gilbert
(LLG) equation yields a Jacobian matrix `B` that maps tangent perturbations
`\delta\bar{m}` to their time derivatives:

```math
\frac{d\bar{m}}{dt} = B \, \bar{m}
```

where `\bar{m}` is the 2N-component tangent vector (two degrees of freedom
per spin site, since `|\mathbf{m}|=1` removes the radial direction). The
eigenvalues `\sigma = \sigma_r + i\sigma_i` of `B` have the following
physical meaning:

- **Imaginary part** `\sigma_i` → oscillation frequency `f = |\sigma_i| / (2\pi)`
- **Real part** `\sigma_r` → damping rate `\Gamma = -\sigma_r / (2\pi)`
  (`\sigma_r < 0` for stable equilibria, so `\Gamma > 0`)

The **least-damped modes** are those with the smallest `|\sigma_r|`, i.e. the
slowest-decaying modes. They dominate the long-time dynamics and are the most
readily observed in experiments. In Arpack, use `which=:LR` (largest real part)
to target them — since `\sigma_r` is negative, "largest real part" means
closest to zero, hence least damped. For the lowest-frequency modes use
`which=:SM` (smallest magnitude), which requires shift-invert via sparse
materialisation.

## Theory

Linearising the Gilbert-form LLG equation

```math
\frac{d\mathbf{m}_i}{dt} = -\gamma \, \mathbf{m}_i \times \mathbf{H}_i^{\rm eff}
+ \alpha \, \mathbf{m}_i \times \frac{d\mathbf{m}_i}{dt}
```

around an equilibrium state `\mathbf{m}_0` (where
`\mathbf{m}_0 \times \mathbf{H}_0 = 0`) gives, for each site `i`:

```math
\delta\dot{\mathbf{m}}_i = -\frac{\gamma}{1+\alpha^2} \Big[
  \mathbf{m}_0 \times \delta\mathbf{H}_i
  + \alpha \, \mathbf{m}_0 \times (\mathbf{m}_0 \times \delta\mathbf{H}_i)
  + \delta\mathbf{m}_i \times \mathbf{H}_0
  + \alpha \, \delta\mathbf{m}_i \times (\mathbf{m}_0 \times \mathbf{H}_0)
\Big]
```

where `\delta\mathbf{H}_i = \sum_j (\partial \mathbf{H}_i / \partial \mathbf{m}_j)\, \delta\mathbf{m}_j`
is the field variation caused by the perturbation. The perturbation
`\delta\mathbf{m}` lives in the local tangent plane perpendicular to
`\mathbf{m}_0`, so we project to a 2D frame per site using rotation matrices
`R_i` and their inverses `R_i^{-1}`.

## Precision Setup

Eigenmode analysis uses symbolic perturbation types (`Epsilon`/`FlatTerm`) for
exact linearisation of local interactions (no finite-difference error). This
requires `AbstractFloat` precision **before** constructing the simulation:

```@example
using MicroMagnetic
using LinearAlgebra
using CairoMakie

MicroMagnetic.set_backend("cpu")
MicroMagnetic.set_precision(AbstractFloat)
nothing #hide
```

## API

### `build_matrix`

```@docs
build_matrix
```

The function supports three output modes:

| Mode        | Keyword           | Return type       | Memory | Use case                     |
| :---------- | :---------------- | :---------------- | :----- | :--------------------------- |
| Dense       | *(default)*       | `Matrix{Float64}` | O(4N²) | Small systems, full spectrum |
| Sparse      | `sparse=true`     | `SparseMatrixCSC` | O(NZ)  | Sparse interactions          |
| Matrix-free | `matrixfree=true` | `LLGJacOperator`  | O(N)   | Large systems, Arpack        |

### `LLGJacOperator`

```@docs
LLGJacOperator
```

### `build_demag_matrix`

```@docs
build_demag_matrix
```

## Example 1: Cubic Anisotropy FMR

We compute the ferromagnetic resonance (FMR) frequency of a single-cell
particle with cubic anisotropy, and compare with the analytical solution.
The setup follows [`test/eigen/test_cubic.jl`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/test/eigen/test_cubic.jl).

For a cubic anisotropy `K_c` and saturation magnetization `M_s`, the
effective anisotropy field along the \[100] direction gives the FMR frequency:

```math
f_{[100]} = \frac{\gamma}{2\pi} (H + K), \quad K = \frac{4 K_c}{\mu_0 M_s}
```

```@example
function setup_cubic(; m0=(0,0,1), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=2e-9)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5)
    add_cubic_anis(sim, 2e4)
    add_zeeman(sim, H)
    init_m0(sim, m0)
    return sim
end

function fmr_frequency_100(H0)
    sim = setup_cubic(m0=(1,0,0), H=(H0, 0, 0))
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
    return imag(eigvals(B)[2]) / 1e9 / (2*pi)
end

function analytical_100(H; Kc=2e4, Ms=8e5)
    gamma = 2.21e5
    K = 4*Kc/(mu_0*Ms)
    return gamma*(H .+ K) / 1e9 / (2*pi)
end

H = 1e5
f_num = fmr_frequency_100(H)
f_ana = analytical_100(H)

println("Numerical:  $(round(f_num, digits=4)) GHz")
println("Analytical: $(round(f_ana, digits=4)) GHz")
println("Relative error: $(round(abs(f_num - f_ana)/f_ana, sigdigits=4))")
```

### Frequency vs. field sweep

By scanning the external field, we reproduce the FMR dispersion shown in
[`figures/eigen/fig2`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/figures/eigen/fig2/fmr.jl):

```@example
Hs = range(0.8e5, 2e5; length=20)

freqs_100 = [fmr_frequency_100(H) for H in Hs]

function fmr_frequency_110(H0)
    sim = setup_cubic(m0=(1,1,0), H=(H0/sqrt(2), H0/sqrt(2), 0))
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
    return imag(eigvals(B)[2]) / 1e9 / (2*pi)
end

freqs_110 = [fmr_frequency_110(H) for H in Hs]

function analytical_110(H; Kc=2e4, Ms=8e5)
    gamma = 2.21e5
    K = Kc/(mu_0*Ms)
    return gamma*sqrt.(H.^2 .- 2K.*H .- 8K^2) / 1e9 / (2*pi)
end

fig = Figure(; size=(500, 360), fontsize=18)
ax = Axis(fig[1, 1]; xlabel="μ₀H (T)", ylabel="Frequency (GHz)")

scatter!(ax, Hs.*mu_0, freqs_100; marker=:rect, markersize=10,
         strokewidth=1, color=:white, label="[100]")
scatter!(ax, Hs.*mu_0, freqs_110; markersize=10,
         strokewidth=1, color=:white, label="[110]")

lines!(ax, Hs.*mu_0, analytical_100.(Hs); color=:slateblue1, label="Analytical")
lines!(ax, Hs.*mu_0, analytical_110.(Hs); color=:slateblue1)

axislegend(ax; position=(0.95, 0.05), labelsize=14)
save("eigen_fig2.png", fig)
nothing #hide
```

![](eigen_fig2.png)

## Example 2: Large-Scale FMR with Matrix-Free Arpack

For large meshes where the dense `2N \times 2N` Jacobian does not fit in
memory, pass `matrixfree=true` to obtain an `LLGJacOperator` and solve with
`Arpack.eigs`. The operator implements the full `AbstractMatrix` interface
needed by Arpack (`size`, `eltype`, 3- and 5-argument `mul!`, `issymmetric`,
`copy`).

This workflow reproduces the standard-problem-4 FMR mode pattern shown in
[`figures/eigen/fig4`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/figures/eigen/fig4/fmr_std.jl):

```julia
using MicroMagnetic
using SparseArrays
using Arpack
using CairoMakie
using Printf

# Step 1: relax to equilibrium using Float64 precision
MicroMagnetic.set_precision(Float64)

function setup_std4(; m0=(0,0,1))
    mesh = FDMesh(nx=24, ny=24, nz=2, dx=5e-9, dy=5e-9, dz=5e-9)
    v = (1, 0.7, 0)
    H = 80e3 .* v ./ sqrt(sum(v.^2))
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5)
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    add_zeeman(sim, H)
    init_m0(sim, m0)
    return sim
end

sim = setup_std4()
relax(sim, stopping_dmdt=0.0001)
m0 = Array(sim.spin)

# Step 2: rebuild with AbstractFloat for eigenmode analysis
MicroMagnetic.set_precision(AbstractFloat)
sim = setup_std4(m0=m0)

# Step 3: build matrix-free operator (O(N) storage)
op = build_matrix(sim, gamma=2.211e5, matrixfree=true)

# Step 4: solve eigenvalues with Arpack
#   which=:SM  → lowest-frequency modes (requires shift-invert via sparse)
#   which=:LR  → least-damped modes (slowest decay, most observable);
#               call eigs(op) directly on the matrix-free operator.
function compute_eigenvalues(op, sim; nev=10, which=:SM)
    if which == :SM
        A = SparseArrays.sparse(op)
        vals, vecs = Arpack.eigs(A; nev=2*nev, which=:SM, tol=1e-6,
                                 maxiter=10*sim.n_total, explicittransform=:none)
    else
        vals, vecs = Arpack.eigs(op; nev=2*nev, which=which, tol=1e-6,
                                 maxiter=10*sim.n_total, explicittransform=:none)
    end
    idx = findall(x -> imag(x) > 0, vals)
    return vals[idx], vecs[:, idx]
end

freqs, vecs = compute_eigenvalues(op, sim, nev=10, which=:SM)

# Step 5: reconstruct 3D magnetization perturbation from tangent eigenvectors
function reconstruct_modes(vecs, sim)
    N = sim.n_total
    mesh = sim.mesh
    nkeep = size(vecs, 2)
    evecs = reshape(vecs, (2, N, nkeep))
    m3d = zeros(ComplexF64, 3, N, nkeep)
    m0 = reshape(sim.spin, 3, N)
    for i in 1:N
        R = MicroMagnetic.rotation_matrix(m0[1,i], m0[2,i], m0[3,i])
        for j in 1:nkeep
            m3d[:, i, j] = R * evecs[:, i, j]
        end
    end
    return reshape(m3d, 3, mesh.nx, mesh.ny, mesh.nz, nkeep)
end

m = reconstruct_modes(vecs, sim)
freqs_GHz = imag.(freqs) ./ (2pi * 1e9)

# Plot |m_y| of the first few modes
fig = Figure(size=(600, 290), fontsize=14)
for i in 1:min(5, length(freqs_GHz))
    vv = abs.(m[2, :, :, 1, i])
    title = @sprintf("%.3f GHz", freqs_GHz[i])
    ax = Axis(fig[1, i]; title=title, aspect=1)
    hidedecorations!(ax)
    heatmap!(ax, vv; colormap=:coolwarm, interpolate=true)
end
save("eigen_fig4.png", fig)
```

!!! note "Matrix-free operator internals"
The `LLGJacOperator` applies the Jacobian in `O(N)` memory using:

```
1. **Unwrap** the 2N tangent vector `x` into a 3N Cartesian `dm` via
   precomputed rotation matrices `R_i`.
2. **Local field variation** `δH_loc`: a single symbolic pass with one
   `Epsilon` tag over `m0 + ε·dm`, extracting the ε-coefficient of the
   effective field. This gives `(∂H_loc/∂m)·dm` in one sweep.
3. **Long-range field variation** `δH_demag`: one `effective_field` call
   per demagnetization interaction (FFT-accelerated).
4. **LLG cross products**: `m0×δH`, `dm×H0`, and damping terms
   `α·m0×(m0×δH)` etc., using cached baseline fields `H0` and `m0×H0`.
5. **Project** the 3D result back to the 2N tangent frame via `R_i^{-1}`.

Compared to the dense builder (which calls `mul!` 2N times to materialise
all columns), the matrix-free `mul!` needs only **1 local symbolic pass +
1 FFT demag call** per application.
```

!!! warning "Restore the precision afterwards"
    `set_precision(AbstractFloat)` changes a global setting. Since `AbstractFloat`
    is not a concrete type, simulations created afterwards will fail (in particular
    on the GPU). When you are done with the eigenmode analysis, restore the default:

```@example
MicroMagnetic.set_precision(Float64)
nothing #hide
```

## References

- X. Fan, S. Zhang, W. Wang, L. Kong, Y. Guo, Y. Liu, and H. Du, "Automatic
  eigenvalue method in micromagnetic and atomistic simulations,"
  *Journal of Magnetism and Magnetic Materials* **622**, 172989 (2025).
  [doi:10.1016/j.jmmm.2025.172989](https://doi.org/10.1016/j.jmmm.2025.172989)

