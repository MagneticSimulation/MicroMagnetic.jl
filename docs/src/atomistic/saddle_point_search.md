```@meta
ShareDefaultModule = true
```

# Saddle-point search (SPS)

Systematic saddle-point search (SPS) starts from one relaxed magnetic state.
It follows low Hessian modes to first-order saddle points, relaxes both sides of
each saddle, and refines accepted connections with the geodesic nudged elastic
band (GNEB) method. Transition mechanisms are identified from the resulting
states and paths; they are not assigned to Hessian modes in advance.

This example uses the atomistic skyrmion model from Müller *et al.* with an
open ``40 × 40 × 1`` lattice, ``J = 1 meV``, ``D/J = 0.45``, and
``B_z = 0.8 D²/(μ_s J)``. The complete runnable script is
[`examples/skyrmion_sps.jl`](../../../examples/skyrmion_sps.jl).

## Define the model

```julia
using MicroMagnetic
using CairoMakie

@using_gpu()

const N = 40
const J = 1.0meV
const D = 0.45J
const MU_S = mu_B
const B_Z = 0.8D^2 / (MU_S * J)

function create_sim()
    mesh = CubicMesh(; nx=N, ny=N, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", save_data=false)
    set_mu_s(sim, MU_S)
    add_exch(sim, J)
    add_dmi(sim, D)
    add_zeeman(sim, (0, 0, B_Z))
    return sim
end
```

Initialize and relax one skyrmion before starting SPS. The example script
contains the short `initial_skyrmion` function used below.

```julia
sim = create_sim()
init_m0(sim, initial_skyrmion)
relax(sim; max_steps=50_000, stopping_dmdt=1e-5,
      save_data_every=-1, save_m_every=-1)
minimum = Float64.(Array(sim.spin))
```

## Search and plot transition paths

```julia
result = find_transitions(
    create_sim, minimum;
    n_modes=8,
    directions=(-1, 1),
    exploration_depth=1,
    n_transitions=3,
    images=21,
    output_folder="skyrmion_sps",
)

plot_transition_paths(
    result; output="skyrmion_sps/transition_paths.png")
```

`find_transitions` records failed candidates and continues exploring the
remaining mode/direction pairs. An accepted path has two distinct endpoint
minima, an index-one saddle, and converged ordinary and climbing-image GNEB
bands. Energies and forces use joules internally; the path plot displays
relative energy in meV by default.

The root output folder contains `minimum.ovf`, `hessian_modes.txt`,
`attempts.txt`, and `transitions_energy.txt`. Each accepted `transition_XX`
folder contains the saddle, 21 path-image OVF files, `energy.txt`, and
`distance.txt`; rejected candidates do not create folders. As in the legacy
`NEB` output, the energy table stores absolute energies in joules and the
distance table stores dimensionless geodesic distances. The combined table
additionally provides energy above the path minimum in joules and meV, while
plots use ``(E-E_min)/meV``, matching the current NEB documentation example.
Set `save_diagnostics=true` only when MMF and GNEB iteration histories are
needed. A CUDA, AMDGPU, oneAPI, or Metal package can be loaded before
`@using_gpu()`; otherwise the same script runs on the CPU.

!!! note
    The full ``40 × 40`` search is a production example and is not executed
    during every documentation build. Its numerical result depends on the
    physical model and tolerances, but users normally only need to change the
    model parameters and the arguments of `find_transitions`, not the solver
    strategy.

## References

- G. P. Müller, P. F. Bessarab, S. M. Vlasov, F. Lux, N. S. Kiselev,
  S. Blügel, V. M. Uzdin, and H. Jónsson, "Duplication, Collapse, and Escape
  of Magnetic Skyrmions Revealed Using a Systematic Saddle Point Search
  Method," *Physical Review Letters* **121**, 197202 (2018).
  [doi:10.1103/PhysRevLett.121.197202](https://doi.org/10.1103/PhysRevLett.121.197202)

- P. F. Bessarab, V. M. Uzdin, and H. Jónsson, "Method for finding mechanism
  and activation energy of magnetic transitions, applied to skyrmion and
  antivortex annihilation," *Computer Physics Communications* **196**,
  335-347 (2015).
  [doi:10.1016/j.cpc.2015.07.001](https://doi.org/10.1016/j.cpc.2015.07.001)

- G. Henkelman and H. Jónsson, "Improved tangent estimate in the nudged
  elastic band method for finding minimum energy paths and saddle points,"
  *Journal of Chemical Physics* **113**, 9978-9985 (2000).
  [doi:10.1063/1.1323224](https://doi.org/10.1063/1.1323224)

- G. P. Müller *et al.*, "Spirit: Multifunctional framework for atomistic
  spin simulations," *Physical Review B* **99**, 224414 (2019).
  [doi:10.1103/PhysRevB.99.224414](https://doi.org/10.1103/PhysRevB.99.224414)

- H. Schrautzer, M. Sallermann, P. F. Bessarab, and H. Jónsson,
  "Identification of mechanisms of magnetic transitions using an efficient
  method for converging on first-order saddle points," *Physical Review B*
  **112**, 104433 (2025).
  [doi:10.1103/z673-hhnp](https://doi.org/10.1103/z673-hhnp)
