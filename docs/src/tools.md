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

## Lorentz TEM (LTEM)

The LTEM tools simulate what a Lorentz transmission electron microscope
would see for a given magnetization configuration: the magnetic phase shift
and, optionally, the Fresnel defocus image.  They also cover **tilted
sample** geometries.  References: Walton *et al.*, J. Phys. D **46**, 175005
(2013) (MALTS); Keimpema *et al.* (2006) (mean inner potential phase);
tilt/projection method after [MagRecon/maglab](https://github.com/MagRecon/maglab).

A worked example that simulates the Fresnel images of a magnetic vortex and
of Bloch/Néel skyrmions is given in
[LTEM Imaging of Typical Magnetic Structures](ltem.md).

### Magnetic phase

For a beam along `+z` the magnetic phase depends only on the beam-integrated
in-plane magnetization `T = ∫ M⊥ dz` and is evaluated as a linearly padded
FFT convolution (exact for the sampled kernel, cross-validated against a
direct real-space sum):

```julia
ovf = read_ovf("skx.ovf")
phi = compute_magnetic_phase(ovf; Ms=8e5)                     # beam along z
phi = compute_magnetic_phase(ovf; Ms=8e5, ty=deg2rad(20))     # tilt about y
phi = compute_magnetic_phase(m, deg2rad(20), "y"; Ms=8e5)     # array method
```

`Ms` is in A/m (the default `1/mu_0` gives the phase of a unit induction),
the returned `N × N` phase is in rad and centered on the sample.  Tilts are
Euler angles in radians about the `x`/`y`/`z` axes (`tx`, `ty`, `tz`); the
volume — magnetization vectors included — is rotated by sequential bilinear
slice rotations while the beam stays fixed.  If the data uses another beam
axis, declare it with `axis` (`:x`, `:y` or `:z`, default `:z`) and the
volume is rotated into place before the tilts are applied.  `N` overrides
the size of the cubic rotation/output grid (by default a power of two that
contains the rotated bounding box).

!!! note "Tilt convention"
    The projected in-plane components are fed to the phase kernel in the
    *laboratory* frame.  Rotating them back into the sample frame (as
    MagRecon/maglab does) mixes in the beam-parallel component: a uniformly
    in-plane magnetized film would get ``1/\cos\alpha`` too much contrast
    and a perpendicularly magnetized film exactly zero contrast at any tilt.
    The laboratory-frame convention reproduces the exact Lorentz deflection
    ``\int \mathbf{B}_\perp\,\mathrm{d}z`` — tilt-independent for the
    in-plane film, ``\propto \tan\alpha`` for the perpendicular one.

### Full defocus image

[`LTEM`](@ref) adds the electric (mean inner potential) phase over the
projected material footprint and propagates the exit wave with the Fresnel
defocus transfer function:

```julia
phase, image = LTEM("skx.ovf"; Ms=8e5, V=300, df=200, ty=deg2rad(20))
```

with `V` the accelerating voltage in kV, `df` the defocus in µm, `V0` the
mean inner potential in V and `alpha` the beam divergence semiangle in rad.
`LTEM` returns the unwrapped magnetic phase and the normalized intensity,
both `N × N` and centered on the sample.  At `df=0` the intensity of a pure
phase object is identically one — a convenient self-check.

### Helpers

`rotate3d`, `warp2d`, `euler_matrix`, `project3d`,
`vector_field_projection`, `vector_padding`, `electron_wavelength`,
`interaction_constant` and `compute_electric_phase` are available for
standalone use (unexported); see the LTEM section of `api.md`.

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
