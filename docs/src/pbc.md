# Periodic Boundary Conditions

MicroMagnetic supports periodic boundary conditions (PBC) on any combination of the
three axes. Periodicity is declared **once, on the mesh**, and every interaction
term respects it automatically: exchange and DMI through the periodic neighbor
table, and the demagnetization field through a dedicated periodic solver. This
removes the silent inconsistency of "periodic exchange but open demag".

## Declaring periodicity

```julia
mesh = FDMesh(; nx=64, ny=64, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9, pbc="xy")
```

`pbc` accepts any combination of the characters `"x"`, `"y"`, `"z"` (default
`"open"`). The rest of the setup is unchanged:

```julia
sim = Sim(mesh; driver="LLG", name="skx_lattice")
set_Ms(sim, 8e5)
add_exch(sim, 1.3e-11)
add_dmi(sim, 3e-3; type="interfacial")
add_demag(sim)                 # solver selected automatically from mesh.pbc
```

## Demagnetization solvers

The demag solver is selected automatically from the mesh periodicity:

| mesh periodicity | solver | method | reference |
| :--- | :--- | :--- | :--- |
| open | open demag | padded-FFT (or direct with `fft=false`) | Newell |
| one axis (`x`/`y`/`z`) | pbc1d | image sums + analytic tail + DC-column fix | Lebecki, *J. Phys. D* **41**, 175005 (2008) |
| two axes (`xy`/`xz`/`yz`) | pbc2d | image sums + analytic tail + DC-column fix | Wang et al., *Comp. Mater. Sci.* **49**, 84 (2010) |
| three axes (`xyz`) | pbc3d | spectral projector (tin-foil convention) | Mansuripur (1989) |

The periodic kernel is the image sum `N^PBC(d) = Σ_images N(d + p·L)`. For one and
two periodic axes this sum converges absolutely, so explicit images plus a
closed-form continuous tail give a well-defined approximation. For three periodic
axes the sum is only conditionally convergent; pbc3d uses the standard spectral
projector instead (the uniform mode carries no field).

### Image counts `Ic` / `Jc`

For pbc1d/pbc2d the number of explicit images per periodic axis defaults to
`Ic = Jc = 4`. Increase it for higher accuracy at linear init-cost:

```julia
add_demag(sim; Ic=8, Jc=8)   # more images, smaller tail error
```

Typical field errors against a brute-force image sum: pbc1d `Ic=2` → 6×10⁻⁵,
`Ic=4` → 2.5×10⁻⁶; pbc2d `Ic=2` → 1×10⁻⁴.

### Legacy explicit-macro mode

The legacy truncated-image macro cell mode is still available and takes precedence
when its keywords are given:

```julia
add_demag(sim; macroPBC=true, Nx=2, Ny=2)   # Nx/Ny/Nz = repetition counts
add_demag(sim; Nx=4, Ny=4)                  # legacy semantics without macroPBC
```

!!! note
    `fft=false` degrades to the direct solver with truncated images; true PBC
    requires the FFT pipeline (default).

## Consistency guarantees

- The periodic axes live on the mesh only; `add_exch`, `add_dmi` and `add_demag`
  all read them from `sim.mesh`. There is no per-term periodicity keyword to set.
- Passing an explicit demag keyword that contradicts the mesh periodic axes (for
  example `Nx` on an open axis) triggers a warning. Passing them on an **open**
  mesh is silent — that is the supported way to test the legacy macro solver in
  isolation.
- The index bookkeeping `indexpbc(i, j, k, nx, ny, nz, xperiodic, yperiodic)`
  wraps around periodic axes; use it instead of `index` when scanning neighbors
  of a periodic mesh by hand.

!!! info "Implementation notes"
    The math, sign conventions and verification methodology of the periodic
    demag solvers are documented for developers in
    [Demag Internals](dev/demag.md).
