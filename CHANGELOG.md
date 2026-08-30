# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## Version [v0.6.0] - unreleased

### Added

- `sim_with`: new `stt` and `sot` keywords forward to `add_stt`/`add_sot` (e.g.
  `stt=(model=:zhang_li, b=-72.438, J=(1,0,0), xi=0.05)`); the torque is applied only in
  stages running with an LLG-family driver so SD minimization never sees it

### Changed

- `sim_with` now validates all keywords up front: unknown keys throw an error (with a
  "did you mean" hint) before the simulation starts, and keys that cannot apply to the
  configured task(s) produce an early warning
- `sim_with` no longer modifies the caller's `Dict`/`NamedTuple` (an internal copy is used)
- `sim_with`: sweeping via `_s`/`_sweep` is now restricted to `task`, `driver`, `Ms` and
  `H`; sweeping any other key throws an error instead of silently pinning its first value
- `sim_with`: `driver_s` is now honored for Relax stages as well, and driver parameters
  (`alpha`, `gamma`, `tol`, ...) are re-applied whenever the driver switches between stages
- `relax`/`run_sim`: `save_data_every=0`/`save_m_every=0` now reliably save the final
  state, including when `max_steps` is reached without convergence; the ovf/vts/vtu
  output branching is unified through a shared internal helper
- `create_sim` copies the keyword `Dict` instead of consuming it, warns about material
  keys ignored by the given mesh type (`mu_s`/`J` on FDMesh, `Ms`/`A` on atomistic
  meshes) and about orphaned keys (`dmi_type` without `D`, `axis` without `Ku`), and its
  unsupported-mesh error message is now informative
- `hysteresis` now returns the collected `(H_values, m_values)`; `plot_hysteresis(H, m)`
  is implemented in the CairoMakie extension to visualize them
- the public API no longer exports `dynamic_matrix` (it was declared but never
  implemented)
- exported the type aliases `NumberOrArrayOrFunction`, `NumberOrTupleOrArrayOrFunction`
  and `TupleOrArrayOrFunction`, which appear in public method signatures
- `run_sim` for `MonteCarlo` is renamed to `run_mc` (the old name still works with a
  deprecation warning)
- `add_mel` is deprecated in favor of `add_magnetoelastic`
- clearer sweep errors: scalar values passed to a `_s` keyword now raise an informative
  `ArgumentError`, and length mismatches report `_s`/`_sweep` instead of `_range`

### Fixed

- fix LaTeX formula rendering in the docs: migrate `docs/src/.vitepress/config.mts` to
  the DocumenterVitepress mathjax plugin and pin the npm dependencies with a lockfile
- fix `ovf2movie` default output filename (undefined `path` variable)
- fix docs/tutorials (std4, std5, skyrmion_stt): `dynamic_m_interval` → `dynamic_m_every`
  and `jld2movie` (removed) → `ovf2movie` on the ovf output folder

### Removed

- the `LLG_STT`, `LLG_CPP` and `LLG_STT_CPP` drivers (deprecated since v0.5.0); use the
  `LLG` driver together with `add_stt`/`add_sot` instead
- `set_ux`/`set_uy`/`set_uz` (driver-bound STT helpers); pass `b`/`J` to `add_stt` instead
- `sim_with` keywords `beta`, `ux`, `uy`, `uz` and `ufun`; put them inside the
  `stt=(model=:zhang_li, ...)` tuple instead
- `create_sim` is no longer exported; build simulations with `Sim` + `set_*`/`add_*`
  (or use the high-level `sim_with`). `MicroMagnetic.create_sim(...)` still works and
  remains the factory used by the NEB/transition machinery

## Version [v0.5.0] - 2026-05-30

### Added

- Web based GUI
- PlotPanel for interactive data visualization
- Color support for isosurface visualization
- Hopfion initialization function and example
- More GUI snippets (80+ code templates)

### Changed

- change relax_m_interval to relax_m_every in sim_with function
- change dynamic_m_interval to dynamic_m_every in sim_with function  
- change relax_data_interval to relax_data_every in sim_with function
- change gif to mp4 for video output (ovf2movie)
- fix GUI arrow sampling order (apply options before calculate positions)
- fix NEB LLG driver initialization (set_initial_condition!)
- fix NEB spin sync for AdaptiveRK integrator

### Removed

- add_exch(sim::MicroSim, geo::Shape, A::Number; name="exch")

### Fixed

- NEB acos error (#22)
- hysteresis loop calculation
- arrow direction in visualization
- sampling adjustment for large meshes

---

## Version [v0.4.0]

### Added
- voronoi
- FEM (exch, demag)
- add_stt, add_sot
- inertia llg

### Changed
- she_torque -> sahe_torque
- EnergyMinimization -> SD

### Removed
- JLD2 

---

## Version [v0.3.9] - 2025-02-16

### Added

- Eigen methods (not finished yet)
- Hexagonal anisotropy
- SHE torque
- Biquadratic exchange

### Fixed

- mu0 for stochastic field