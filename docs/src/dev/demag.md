# Demag Internals

```@meta
CurrentModule = MicroMagnetic
```

The demagnetization field is the most expensive term per step and the only one
with a global (long-range) coupling, so it has its own solver family. The user
level of this — how to turn periodicity on — is on the
[Periodic BC page](../pbc.md); this page records the implementation, math and
verification methodology.

## Solver family

| periodicity | solver (file) | method |
| :--- | :--- | :--- |
| open | `micro/demag.jl`, `micro/demag_direct.jl` | padded-FFT convolution with Newell kernels (or direct summation, `fft=false`) |
| one axis | `micro/demag_pbc1d.jl` | image sums + analytic 1D tail + DC-column fix (Lebecki 2008) |
| two axes | `micro/demag_pbc2d.jl` | image sums + analytic 2D tail + DC-column fix (Wang 2010) |
| three axes | `micro/demag_pbc3d.jl` | spectral projector +k̂k̂, tin-foil (Mansuripur 1989) |
| legacy | `micro/demag.jl` | truncated image sums ("macro"), kept for compatibility |

The periodic kernel is `N^PBC(d) = Σ_images N(d + p·L)`. 1D and 2D image sums
converge absolutely (1/r³ decay with dipolar cancellation), so explicit images
plus a closed-form continuous tail is a strict decomposition. The 3D sum is only
conditionally convergent and has no continuous tail — hence the spectral
projector, under which the uniform mode carries no field.

## The open/pbc1d/pbc2d pipeline (macro route)

pbc1d and pbc2d **reuse the entire padded-FFT pipeline** of the open solver
(buffer allocation, FFT plan, packing, spectral multiply, in-place c2r,
collect). Periodicity is implemented as three post-processing steps on the
quadrant kernels:

1. **Image sums** — the tensor-building kernels take non-zero image counts
   (`Ic`, `Jc`);
2. **Per-entry analytic tail** — a KernelAbstractions kernel adds the
   closed-form continuous-tail contribution (1D: Lebecki's continuous
   approximation `T_ab(x0, v, w)`, error ~Xc⁻³; 2D: Wang 2010 Eq. 10–15
   inclusion–exclusion `Ninf`, error ~γ^(−3/2));
3. **DC-column fix** — the uniform-mode response is replaced by its analytic
   value (1D: rod identity, exactly zero for axial uniform states; 2D: zero for
   in-plane uniform states, exact self-term +1 for open-axis diagonal
   components).

Why padded FFT gives periodicity for free: the circular wrap of the padded
convolution implements the image sum exactly on the periodic axes (the kernel is
already periodic), while zero-padded magnetization degrades open axes to the
linear (open) response. The wrap never touches the Nyquist plane, so zeroing it
is harmless.

pbc3d is deliberately separate: the projector is direction-independent,
initialisation is free (no Newell evaluation) and each step is ~2× faster than
the macro route.

## Sign conventions (classic trap)

The package kernel is the **negative** of the field kernel: `H = -Σ N·M`, so
far-field `N_xx(u, 0, 0) < 0`. Self-derived tail formulas written in the field-kernel
convention (positive `(2u²−ρ²)/r⁵` structure) are off by a global sign. The
published `Ninf` (Wang 2010) and `N^cont` (Lebecki 2008) are both already in the
package convention — no negation. Cross-checking magnitudes but not signs has
burned a full debugging cycle here once; compare signs too.

Two more traps from implementation history, both caught by brute-force baselines:

- **Permutation tables**: kernel components are built in a canonical frame and
  mapped back (`PBC1D_PCOMP` / `PBC2D_PCOMP`). A transcription error left axial
  components correct while lateral tensor DC columns got garbage — random
  textures looked fine, so cross-checks must cover every axis/direction pair.
- **Test-state components**: `m[3*(axis-1)+1:3:end]` starts at the wrong offset
  for axis=2 (sets x, not y). The lateral-uniform test state has a *legitimate*
  nonzero response, so a "wrong but physical-looking" signal passes casually.
  Use `m[axis:3:end]`.

## Verification methodology

The pattern that works, in order:

1. **Brute-force image-sum baseline first**: direct summation with the same sign
   convention, explicit images over all periodic directions. Self-check the
   baseline itself: macro-N ≡ brute-force-N at 1e-11 proves the baseline before
   any comparison.
2. **Per-buffer pipeline diff**: `m_pad → M_pad → M_out → h_pad` compared stage
   by stage against a hand-written Julia pipeline; the first divergent buffer
   localises the bug.
3. **Zero-response probes**: states that must produce exactly zero field (1D
   axial uniform, 2D in-plane uniform) are the most sensitive detectors — any
   error shows up in full. Caveat: they are "small difference of large numbers",
   so they do not represent general accuracy, and the test state itself must be
   verified first.
4. **Analytic special cases** from Lebecki App. B.7 and Wang Eq. 13 as
   independent anchors.

## Accuracy and cost (F64/CPU, measured)

| item | value |
| :--- | :--- |
| pbc1d field error vs brute force | Ic=2 → 6×10⁻⁵, Ic=4 → 2.5×10⁻⁶ |
| pbc2d field error | Ic=2 → ~1×10⁻⁴ |
| uniform-state (amplified) error | pbc1d Ic=2 → 2.2×10⁻⁴·Ms, Ic=6 → ≈0 |
| init (kernel construction) | open 0.07 s, pbc1d 0.57 s (DC fix 0.25 s), pbc2d 0.4–2.4 s (grows with layers), pbc3d ≈0 |
| per step | pbc3d ~2× faster than macro; pbc1d/pbc2d (macro route) = open-route speed |

The FFT plan is cached and reused across steps. GPU demag currently relies on
the cuFFT/default-stream sharing assumption; porting to other GPU backends needs
that assumption re-checked (see [Performance](performance.md)).

The neighbor-loop terms (exchange, DMI) have their own acceleration layer —
see [Parameter System](params.md) for the class/pair-table stencil partition.
