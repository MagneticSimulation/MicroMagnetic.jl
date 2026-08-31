# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## Version [v0.6.0] - unreleased

### Added

- `sim_with`: new `stt` and `sot` keywords forward to `add_stt`/`add_sot` (e.g.
  `stt=(model=:zhang_li, b=-72.438, J=(1,0,0), xi=0.05)`); the torque is applied only in
  stages running with an LLG-family driver so SD minimization never sees it
- new `Fill{T,N}` array type with O(1) storage for uniformly valued parameters (vendored
  implementation, no FillArrays.jl dependency): uniform parameters are no longer materialised
  into O(n) arrays, and an isbits `Fill` is passed to GPU kernels by value, where the per-spin
  read folds into a compile-time immediate. Design principle: "spatial or not" is a property of
  the value, not of the type or the API
- `set_alpha(sim, alpha)` now accepts a plain number in addition to arrays and functions
  (`NumberOrArrayOrFunction`); a uniform alpha is stored as a zero-allocation `Fill`
- new tests `test/test_fill.jl` and `test/test_fill_equivalence.jl`: parameters stored as `Fill`
  and the same values materialised into dense arrays give bit-identical trajectories (`alpha`,
  `Ms`, `Ku` and `D`) on the CPU and CUDA backends
- magnetoelastic coupling via `add_magnetoelastic(sim; ...)`: `model=:tensor` couples the
  magnetisation to a stress tensor (`lambda_s` + `sigma`) and `model=:cubic` uses the cubic
  magnetoelastic constants `B1`/`B2` with a six-component strain `[εxx, εyy, εzz, εxy, εxz, εyz]`;
  `sigma`/`strain` accept a 4D array or an `(i, j, k) -> (...)` function
- twin monoclinic anisotropy interaction `add_twin_mono_anis` (with the new example
  `examples/twin_monoclinic_anisotropy.jl`)
- systematic saddle-point search and transition-path workflow: `find_saddle`,
  `find_transitions`, `compute_hessian_modes` (HessianModes) and GNEB with climbing-image
  support (`relax_gneb!`, `set_climbing_image!`, `set_image_type!`, `gneb_images`);
  `plot_transition_paths` visualises the result in the CairoMakie extension. New docs page
  `atomistic/saddle_point_search.md` and example `examples/skyrmion_sps.jl`
- matrix-free mode for the eigenvalue method: `build_matrix(sim; matrixfree=true)` returns an
  `LLGJacOperator` (O(N) memory; a mat-vec is one symmetry pass plus one demag FFT) that can be
  fed to iterative eigensolvers such as `Arpack.eigs`; the assembled-matrix path is optimised as
  well. New docs page `eigen/eigenmodes.md`
- the eigen matrix assembly accepts symbolic parameters (Symbolics `FlatTerm`) for DMI, with a
  `cpart` accessor for symbolic term components

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
- `create_sim` is deprecated and will be removed in v0.7.0: build simulations with
  `Sim` + `set_*/add_*` (low level, see basics.md) or use the high-level `sim_with`;
  it remains exported and fully functional in v0.6.0
- clearer sweep errors: scalar values passed to a `_s` keyword now raise an informative
  `ArgumentError`, and length mismatches report `_s`/`_sweep` instead of `_range`
- the `SpatialLLG` driver is merged into `LLG`: `LLG.alpha` is now an `AbstractArray{T,1}`
  (by default a zero-allocation `Fill(0.1, n)`), so uniform and spatial damping run through the
  same kernel and call-back; `driver="SpatialLLG"` still works but is deprecated (a deprecation
  warning is emitted and it is treated as `"LLG"`, so output files are named `*_llg.txt` either
  way) and will be removed in v0.7
- all per-spin term parameters (`SpatialExchange.A`, `InterfacialDMI.D`,
  `SpatialVectorBulkDMI.Dx/Dy/Dz`, `Anisotropy.Ku`, `CubicAnisotropy.Kc`,
  `StochasticField.temperature`, `ZhangLiTorque.xi`, `SlonczewskiTorque.J`, `DFTorqueField.aj`,
  ...) are routed through a single internal helper (`make_param`): a uniform number is stored as
  a `Fill`, a function or array is materialised into a dense array
- `sim.mu0_Ms`/`sim.mu_s`/`sim.pins` now default to O(1) `Fill`s (0/false), and
  `set_Ms`/`set_mu_s`/`set_pinning` replace the storage instead of writing into it; `set_Ms`
  keeps the A/m → Tesla (`mu_0`) scaling semantics
- GPU kernels no longer mark parameters that may be a `Fill` (`mu0_Ms`, `Ms`/`mu_s`, `Ku`, `Kc`,
  `A`, `D`, `pins`, ...) with `@Const`, and the stochastic-field kernel reads alpha per spin, so
  a spatially varying alpha now also applies to thermal noise
- demag memory/speed overhaul: the Fourier tensors are packed into the (y,z) parity
  fundamental domain (~4x less memory), the inverse transform uses an in-place c2r plan
  (the padded `h_pad` copy is gone) and the 1/N scaling is folded into the tensors
- `sim.driver.alpha = <scalar>` keeps working with the per-cell damping array: a scalar
  is stored as a uniform O(1) `Fill`, arrays are accepted with the eltype matched to the
  driver precision (length is validated), and function-form parameters go through
  `set_alpha(sim, ...)` (which needs the mesh for spatial sampling)

### Fixed

- fix LaTeX formula rendering in the docs: migrate `docs/src/.vitepress/config.mts` to
  the DocumenterVitepress mathjax plugin and pin the npm dependencies with a lockfile
- fix `ovf2movie` default output filename (undefined `path` variable)
- fix docs/tutorials (std4, std5, skyrmion_stt): `dynamic_m_interval` → `dynamic_m_every`
  and `jld2movie` (removed) → `ovf2movie` on the ovf output folder
- `relax` now works with the `InertialLLG` driver (previously it fell into the minimization
  branch): both `relax(sim::AbstractSim)` and `relax(sim::NEB)` decide between time integration
  and minimization by checking whether the driver carries an integrator instead of testing
  `isa(driver, LLG)`
- fix the interlayer exchange energy calculation
- fix the interlayer DMI energy calculation
- fix symbolic type truncation in the local `H_eff` kernels

### Removed

- the `LLG_STT`, `LLG_CPP` and `LLG_STT_CPP` drivers (deprecated since v0.5.0); use the
  `LLG` driver together with `add_stt`/`add_sot` instead
- `set_ux`/`set_uy`/`set_uz` (driver-bound STT helpers); pass `b`/`J` to `add_stt` instead
- `sim_with` keywords `beta`, `ux`, `uy`, `uz` and `ufun`; put them inside the
  `stt=(model=:zhang_li, ...)` tuple instead
- in-place writes into `sim.mu0_Ms`/`sim.mu_s`/`sim.pins` (`copyto!` or broadcast assignment)
  now throw while these fields hold a `Fill`; use `set_Ms`/`set_mu_s`/`set_pinning`, which
  replace the storage
- the Enzyme dependency (and the `EnzymeExt` package extension) is removed; automatic
  differentiation is no longer used anywhere in the package

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