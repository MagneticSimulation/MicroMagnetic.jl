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

### Fixed

- fix `ovf2movie` default output filename (undefined `path` variable)
- fix docs/tutorials (std4, std5, skyrmion_stt): `dynamic_m_interval` → `dynamic_m_every`
  and `jld2movie` (removed) → `ovf2movie` on the ovf output folder

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