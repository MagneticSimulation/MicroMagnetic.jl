# Lorentz TEM (LTEM) tooling

MicroMagnetic ships a small set of tools to simulate Lorentz transmission
electron microscopy (LTEM) observables — the magnetic phase shift and
Fresnel defocus images — directly from a micromagnetic magnetization
configuration, including **tilted sample** geometries.

```@contents
Pages = ["ltem.md"]
Depth = 2
```

## Physics

For an electron beam propagating along `+z` through a thin sample, the
magnetic contribution to the phase is (electron charge ``q=-e``)

```math
\phi(\mathbf{r}) = \frac{q}{\hbar}\int A_z(\mathbf{r},z)\,\mathrm{d}z,
```

which in the image plane depends only on the beam-integrated in-plane
magnetization ``\mathbf{T}(\mathbf{r}) = \int \mathbf{M}_\perp\,\mathrm{d}z``:

```math
\phi_{\mathbf{k}} = \frac{i\mu_0}{2\Phi_0}\,
\frac{f_y \hat{T}_x - f_x \hat{T}_y}{f^2},
\qquad \Phi_0 = \frac{h}{2e}.
```

The phase is computed by a linearly padded FFT convolution with the sampled
real-space kernel, which reproduces the exact direct real-space sum
(`compute_magnetic_phase_direct`) and therefore the continuum integral.

Tilted samples are handled by **actively rotating the volume** —
magnetization vectors included — while the beam stays fixed along `z`
(sequential bilinear 2-D slice rotations, the method used by
[MagRecon/maglab](https://github.com/MagRecon/maglab)).  The projected
in-plane components are then fed to the kernel in the laboratory frame,
i.e. *without* rotating them back into the sample frame; this reproduces
the exact Lorentz deflection ``\int \mathbf{B}_\perp\,\mathrm{d}z`` (e.g. a
uniformly in-plane magnetized film keeps its contrast under tilt, while a
perpendicularly magnetized film gains contrast ``\propto \tan\alpha``).

The electric (mean inner potential) phase follows
``\phi_E = \sigma V_0 t`` with the interaction constant
``\sigma \approx 6.53\times 10^6\ \mathrm{rad/(V\,m)}`` at 300 kV, and the
Fresnel defocus image is formed with the transfer functions of MALTS
(Walton et al., 2013).

## Magnetic phase

```@example ltem
using MicroMagnetic

# a vortex state on a 24x24x4 film (4 nm cells)
nx, ny, nz = 24, 24, 4
m = zeros(3, nx, ny, nz)
for i in 1:nx, j in 1:ny
    phi = atan(j - (ny+1)/2, i - (nx+1)/2)
    env = exp(-(((i-(nx+1)/2)^2 + (j-(ny+1)/2)^2)/25.0)^4)
    m[1, i, j, :] .= -sin(phi) * env
    m[2, i, j, :] .=  cos(phi) * env
end

phi = compute_magnetic_phase(m; Ms=8e5, dx=4e-9, dy=4e-9, dz=4e-9)
println(size(phi))   # 32x32, centered on the sample
```

An `OVF2` file can be passed directly (voxel sizes are taken from the file):

```julia
ovf = read_ovf("magnetization.ovf")
phi = compute_magnetic_phase(ovf; Ms=8e5)
```

### Tilting the sample

Tilts are given in radians about the `x`, `y` or `z` axes:

```julia
phi  = compute_magnetic_phase(ovf; Ms=8e5, ty=deg2rad(20))                # tilt about y
phi  = compute_magnetic_phase(m, deg2rad(20), "y"; Ms=8e5)                # equivalent
phi  = compute_magnetic_phase(ovf; Ms=8e5, tx=0.1, ty=0.2, tz=deg2rad(5)) # full Euler
```

If the magnetization data is stored with a different beam axis (e.g. the
film normal along `x`), declare it with `axis` and it is rotated into place
before the tilts are applied:

```julia
phi = compute_magnetic_phase(ovf; Ms=8e5, axis=:x)   # beam along the old +x
```

## Full LTEM image

`LTEM` adds the electric phase and propagates the exit wave with the Fresnel
defocus transfer function:

```@example ltem
phase, image = LTEM(m; Ms=8e5, dx=4e-9, dy=4e-9, dz=4e-9,
                    V=300, df=200, V0=0, ty=deg2rad(20))
println(extrema(image))
```

with `V` the accelerating voltage (kV), `df` the defocus (µm), `V0` the mean
inner potential (V) and the beam divergence semiangle `alpha` (rad, default
``10^{-5}``).  At `df = 0` the returned intensity is identically one for a
pure phase object, which makes a convenient self-check.  `LTEM` also accepts
an `OVF2` object or a file name.

## Helpers

The rotation/projection machinery is available on its own:

- `rotate3d(vol, tx, ty, tz)` — rotate a 3-D scalar or `(3, N, N, N)` vector
  field by Euler angles (sequential bilinear slice rotations);
- `warp2d(img, theta)` — rotate a square 2-D array;
- `project3d(v)` — integrate along the beam direction;
- `vector_field_projection(m, theta, axis)` — rotate (degrees) and project;
- `vector_padding(m, N)`, `pad_array(v, dims...)` — centered padding;
- `electron_wavelength(V)`, `interaction_constant(V)`,
  `compute_electric_phase(V, V0, Lz, beta)` — electron optics.

## References

- S. K. Walton *et al.*, "MALTS: A tool to simulate Lorentz Transmission
  Electron Microscopy from micromagnetic simulations",
  J. Phys. D: Appl. Phys. **46**, 175005 (2013).
- L. Keimpema *et al.*, "Electron Holography Image Simulation of
  Nanoparticles" (2006) — mean inner potential phase.
- [MagRecon/maglab](https://github.com/MagRecon/maglab) — tilt/projection
  method (note: maglab additionally rotates the projected vector back into
  the sample frame, which cancels the tilt dependence of the integrated
  induction; the implementation here deliberately does not).
