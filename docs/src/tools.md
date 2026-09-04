# Analysis Tools

MicroMagnetic ships a few analysis routines that operate on simulation objects,
magnetization arrays or OVF files. They live in `src/tools/`.

## Topological charge

[`compute_skyrmion_number`](@ref) computes the skyrmion number *Q* (topological
charge) with the Berg–Lüscher method. It works on a flat magnetization array
(`sim.spin` is exactly that) or directly on an OVF file:

```julia
Q = compute_skyrmion_number(Array(sim.spin), mesh)   # current simulation state
Q = compute_skyrmion_number("snapshot.ovf")          # any OVF file on disk
```

The per-layer variant `compute_skyrmion_number_layers("snapshot.ovf")` returns a
vector with the charge of each layer — useful for films where the skyrmion number
per slice matters.

## Guiding center

[`compute_guiding_center`](@ref) returns the guiding center `(Rx, Ry)` of a
vortex or skyrmion (Thiele coordinate). It accepts a simulation or a flat
magnetization array, with optional windowing over a sub-box:

```julia
Rx, Ry = compute_guiding_center(sim)
Rx, Ry = compute_guiding_center(sim; z=3)   # restrict to layer k=3
```

To track the guiding center over a whole run, add it as an extra column of the
saver instead of post-processing snapshots:

```julia
item = SaverItem(("Rx", "Ry"), ("<m>", "<m>"), compute_guiding_center)
push!(sim.saver.items, item)
```

This is how the vortex-core dynamics example (`micromagnetics/skyrmion_stt.md`)
records trajectories.

## Lorentz TEM (LTEM) phase retrieval

[`compute_magnetic_phase`](@ref) computes the magnetic phase shift that an LTEM
would measure for a given magnetization configuration, using a Fourier-space
projection along a chosen axis (Beleggia/Phatak-type approach). The phase is
returned as an `Nx × Ny` matrix in rad:

```julia
ovf = read_ovf("skx.ovf")
phi = compute_magnetic_phase(ovf; Ms=8e5, axis=:z)   # OVF2 container
```

There is also an array method, `compute_magnetic_phase(m, theta, axis; N1, N2,
Ms, d)` for a 4D magnetization array of shape `(3, nx, ny, nz)`, where the
sample can be tilted by the Euler angle `theta` before projection (rotation
axis `"X"` or `"Y"`) and the padding of the projection / Fourier kernel is
controlled with `N1`/`N2`.

!!! note
    The LTEM routines assume a sample of uniform thickness `Lz` taken from the
    mesh step size, and the phase convention follows the under-focus intensity
    transfer commonly used in Lorentz microscopy.

## Voronoi tessellation (polycrystalline samples)

[`voronoi`](@ref) generates a Voronoi tessellation of the mesh and assigns a
region id to every cell — the starting point for polycrystalline simulations
where grain boundaries carry different material parameters:

```julia
mesh = FDMesh(; nx=200, ny=200, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9)
voronoi(mesh; min_dist=20, seed=123456)   # fills mesh.region with grain ids
plot_voronoi(mesh)                        # inspect the tessellation

# per-grain parameters via region_map
set_Ms(sim, region_map(1 => 8.6e5, 2 => 7.9e5))
add_exch(sim, region_map(1 => 1.3e-11, 2 => 1.0e-11))
```

`min_dist` (in cells) controls the minimum grain-center separation, `seed` makes
the tessellation reproducible. Regions created this way compose with
[`set_region`](@ref) and the region-based parameter functions described on the
[basics page](basics.md).
