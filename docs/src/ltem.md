```@meta
ShareDefaultModule = true
```

# LTEM Imaging of Typical Magnetic Structures

This example simulates Lorentz TEM images of three canonical magnetic
structures — a **magnetic vortex**, and a **skyrmion with Bloch or Néel
helicity** — and reproduces the signatures known from experiment and from
LTEM simulation studies:

* the **circular magnetic-step phase** of a vortex and its bright/dark ring
  contrast at defocus;
* the **helicity dependence** of the skyrmion contrast: a Bloch skyrmion is
  visible at zero tilt, while a pure Néel skyrmion is *invisible* at zero
  tilt and only appears once the sample is tilted;
* the **contrast reversal** between underfocus and overfocus.

Import the modules and define two small helpers for plotting:

````@example ltem
using MicroMagnetic
using CairoMakie

function show_phase!(pos, phi, title)
    ax = Axis(pos; title=title, aspect=DataAspect())
    heatmap!(ax, phi; colormap=:balance)
    hidedecorations!(ax)
end

function show_image!(pos, img, title)
    ax = Axis(pos; title=title, aspect=DataAspect())
    heatmap!(ax, img; colormap=:grays, colorrange=(0.85, 1.15))
    hidedecorations!(ax)
end
````

Throughout this page we use a 300 kV microscope and skip the electric
(mean inner potential) phase by setting `V0 = 0`, so that all contrast is
purely magnetic. Defocus values are in µm and tilts in rad.

## Magnetic vortex in a disk

A vortex in a circular disk has an in-plane magnetization that curls around
the center and an out-of-plane core. We build the standard vortex ansatz on
a 320 nm × 320 nm window (cell size 2.5 nm) for a disk of radius 95 nm and
thickness 24 nm:

````@example ltem
n = 128; dx = dy = 2.5e-9; dz = 24e-9
R = 95e-9
m = zeros(3, n, n, 1)
for i in 1:n, j in 1:n
    x = (i - 0.5) * dx - 64e-9
    y = (j - 0.5) * dy - 64e-9
    r = hypot(x, y)
    if r < R
        f = tanh(r / 20e-9)                      # in-plane curl amplitude
        m[1, i, j, 1] = -y / r * f
        m[2, i, j, 1] =  x / r * f
        m[3, i, j, 1] = -tanh(r / 8e-9)          # down core
    end
end
nothing #hide
````

The beam-integrated in-plane induction of a vortex is azimuthal, and the
magnetic phase is therefore a **circular step**: the phase is constant on
each side of the core and drops through the core region (by about −2.5 rad
for this thickness and magnetization). Compute the phase and three Fresnel
images:

````@example ltem
phi_v = compute_magnetic_phase(m; Ms=8e5, dx=dx, dy=dy, dz=dz)
im_v1 = LTEM(m; Ms=8e5, V=300, df=-300, V0=0, dx=dx, dy=dy, dz=dz)[2]
im_v2 = LTEM(m; Ms=8e5, V=300, df=+300, V0=0, dx=dx, dy=dy, dz=dz)[2]
im_v3 = LTEM(m; Ms=8e5, V=300, df=-300, V0=0, ty=0.30, dx=dx, dy=dy, dz=dz)[2]
nothing #hide
````

````@example ltem
fig = Figure(size=(1050, 760))
show_phase!(fig[1, 1], phi_v, "magnetic phase")
show_image!(fig[1, 2], im_v1, "df = -300 um, tilt 0")
show_image!(fig[2, 1], im_v2, "df = +300 um, tilt 0")
show_image!(fig[2, 2], im_v3, "df = -300 um, tilt 17 deg")
fig
````

What we see:

* the **phase map** shows the circular magnetic step centered on the vortex
  core (about −2.5 rad total for this thickness);
* at zero tilt the defocused images show **concentric bright/dark rings**
  around the core position — the Fresnel image of the step — with the
  contrast **reversed between underfocus and overfocus**;
* tilting the sample mixes the out-of-plane core into the beam-integrated
  induction and adds the familiar **asymmetric core contrast**.

## Skyrmions: Bloch versus Néel helicity

For a skyrmion the in-plane induction winds around the core with a fixed
helicity angle ``\chi``: ``\mathbf{m}_\perp \propto \cos\chi\,\hat{r} +
\sin\chi\,\hat{\phi}``. The LTEM phase kernel is blind to the curl-free
(radial) part of the projected induction, so **the zero-tilt contrast of a
skyrmion scales with ``\sin\chi``**: a Bloch skyrmion (``\chi = 90°``)
shows the textbook contrast at zero tilt, while a pure Néel skyrmion
(``\chi = 0°``) has *identically zero* phase at zero tilt and only becomes
visible when the sample is tilted.

We use the built-in [`skyrmion`](@ref) texture (radius 30 nm, wall width
15 nm) on a 240 nm window with 2 nm cells and a film thickness of 50 nm:

````@example ltem
n2 = 120; dx2 = dy2 = 2e-9; dz2 = 50e-9
function skyrmion_m(helicity)
    m = zeros(3, n2, n2, 1)
    sk = skyrmion(; center=(120e-9, 120e-9), R=30e-9, p=-1, v=1, phi=helicity)
    for i in 1:n2, j in 1:n2
        m[:, i, j, 1] .= sk((i - 0.5) * dx2, (j - 0.5) * dy2, 0.0)
    end
    return m
end
m_bloch = skyrmion_m(pi / 2)   # Bloch: in-plane component azimuthal
m_neel  = skyrmion_m(0.0)      # Néel: in-plane component radial
nothing #hide
````

````@example ltem
phi_b = compute_magnetic_phase(m_bloch; Ms=8e5, dx=dx2, dy=dy2, dz=dz2)
b_m = LTEM(m_bloch; Ms=8e5, V=300, df=-300, V0=0, dx=dx2, dy=dy2, dz=dz2)[2]
b_p = LTEM(m_bloch; Ms=8e5, V=300, df=+300, V0=0, dx=dx2, dy=dy2, dz=dz2)[2]
n_0 = LTEM(m_neel;  Ms=8e5, V=300, df=-300, V0=0, dx=dx2, dy=dy2, dz=dz2)[2]
n_t = LTEM(m_neel;  Ms=8e5, V=300, df=-300, V0=0, ty=0.30, dx=dx2, dy=dy2, dz=dz2)[2]
nothing #hide
````

````@example ltem
fig = Figure(size=(1250, 760))
show_phase!(fig[1, 1], phi_b, "Bloch: phase (tilt 0)")
show_image!(fig[1, 2], b_m, "Bloch: df = -300 um, tilt 0")
show_image!(fig[1, 3], b_p, "Bloch: df = +300 um, tilt 0")
show_phase!(fig[2, 1], phi_b .* 0, "Neel: phase (tilt 0) = 0")
show_image!(fig[2, 2], n_0, "Neel: df = -300 um, tilt 0")
show_image!(fig[2, 3], n_t, "Neel: df = -300 um, tilt 17 deg")
fig
````

Reading the figures:

* **Bloch skyrmion (top row)**: the phase has a two-lobed dipolar
  structure, and the Fresnel images show a bright/dark lobe pair that
  **exchanges lobes between underfocus and overfocus**;
* **Néel skyrmion (bottom row)**: the zero-tilt phase and image are
  *identically flat* — nothing would be seen in an untilted microscope —
  while a 17° tilt produces a clear contrast ring around the skyrmion.

This helicity dependence is a well-known result of Lorentz microscopy and
is reproduced here without any adjustable parameter.

## Notes

* The contrast at zero tilt comes exclusively from the in-plane induction;
  a perpendicular magnetization component becomes visible only through
  tilting (or through its 1/r tail at the sample edges).
* All images are normalized (the intensity of the incident beam is 1), so
  bright means *more* electrons than the background and dark means fewer.
* For real samples remember the electric phase (`V0`) adds a strong but
  featureless background over the material footprint.

## References

* J. N. Chapman and M. R. Scheinfein, "Transmission electron microscopies
  of magnetic microstructures", J. Magn. Magn. Mater. **200**, 729 (1999).
* C. Phatak, A. K. Petford-Long and M. De Graef, "Recent advances in
  Lorentz microscopy", Curr. Opin. Solid State Mater. Sci. **20**, 107
  (2016).
* S. K. Walton *et al.*, "MALTS: A tool to simulate Lorentz Transmission
  Electron Microscopy from micromagnetic simulations", J. Phys. D **46**,
  175005 (2013).
