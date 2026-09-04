# Driver Internals

```@meta
CurrentModule = MicroMagnetic
```

The driver is the equation of motion attached to a simulation: energy
minimisation (`SD`), Gilbert dynamics (`LLG`, also used for stochastic LLG via
`add_thermal_noise`), and inertial dynamics (`InertialLLG`). `MonteCarlo` is a
separate `AbstractSim` subtype with its own sweep loop, not a driver.

## Types

```julia
SD{T}          # steepest descent: preconditioned energy minimisation
LLG{T}         # alpha::AbstractArray{T,1} (Fill or dense), gamma::T, integrator
InertialLLG{T} # scalar alpha, gamma, tau (angular-momentum relaxation time)
```

All drivers are parameterised by the element type `T` frozen at `Sim`
construction (see [Architecture](architecture.md#global-state-precision-backend-group-size)).

**`SpatialLLG` no longer exists.** When `alpha` became a (possibly `Fill`) array,
`LLG` absorbed the spatial variant; the string `"SpatialLLG"` is accepted by
`set_driver` only as a deprecation alias that warns and maps to `"LLG"`.

## set_driver: dispatch and guards

`set_driver(sim; driver=..., kwargs...)` is the only supported way to switch a
driver after construction (assigning `sim.driver` directly bypasses the guards).
The flow:

1. **Name normalisation** (`_normalize_driver_name`): `"SpatialLLG"` → `"LLG"`
   with a `depwarn`; unknown names throw with the supported list
   (`"None"`, `"SD"`, `"LLG"`, `"InertialLLG"`).
2. **Type dispatch**: `create_driver(::Type{T}, name, integrator, n_total)`
   constructs the concrete driver struct with the simulation's precision.
3. **Guards** (the A1/A2/A3 rules from the 2026-09 redesign):
   - keyword `alpha` accepts numbers and arrays; a *function* is rejected —
     spatial sampling needs the mesh, use `set_alpha(sim, f)`;
   - the integration tolerance lives on `driver.integrator.tol`
     (`set_driver(sim; tol=...)`); the driver-level `tol` field was removed;
   - `gamma` and other scalar fields are converted to the driver's `T` at
     assignment.

### Property assignment on LLG

`sim.driver.alpha = x` is intercepted by a `setproperty!` hook:

- `x isa Function` → `ArgumentError` (use `set_alpha`);
- `x isa Number` → in-place `Fill(T(x), n)`;
- `x isa AbstractArray` → length must match, eltype converted to `T`.

Other (typed scalar) fields fall through to the default convert-on-assign
semantics, so `driver.gamma = 2.21e5` keeps working in Float32 sessions.

## Integrators

`LLG` and `InertialLLG` integrate through the `Integrator` hierarchy in
`src/integrator/`: adaptive explicit RK (`DOPRI54`, CashKarp, Fehlberg, ...),
Cayley-transformed integrators preserving `|m| = 1` (`RK_Cayley`,
`DOPRI5_Cayley`), Heun and GPSM. Two contracts worth knowing:

- integrator state lives inside the driver; `relax` currently maintains its own
  `prespin` copies around the integrator (a known wart, see
  [Performance](performance.md));
- every RHS evaluation is the full effective-field pipeline, so the RHS count
  per step is the primary cost lever (a FSAL upgrade for DOPRI54 — 7 → 6
  evaluations — is on the roadmap).

## Stopping criteria

`relax` measures `dm/dt` (or the torque norm for `SD`) each step and stops below
`stopping_dmdt` (default 0.01). `run_until(sim, t)` integrates to an absolute
time. Both write through the saver, so a relaxation trace is always inspectable
as a table.

## Adding a driver

Follow the `InertialLLG` precedent: a `Driver` subtype with `T` parameters, a
branch in `create_driver` (plus its RHS kernel and a `call_back` that runs the
effective-field pipeline first), a `set_driver` name entry, and tests
(`test/test_set_driver.jl` shows the guard/regression pattern, including the
`@test_logs` match-mode pitfalls). The Cayley files in `src/llg/` show how a
structure-preserving scheme plugs in at the integrator level instead.
