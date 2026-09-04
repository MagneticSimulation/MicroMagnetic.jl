# Performance Notes

```@meta
CurrentModule = MicroMagnetic
```

Working notes for making things fast without breaking things. The full
issue-by-issue tracker lives in the repository root (`PERFORMANCE.md`); this page
holds the stable knowledge.

## Measure first

- Timing is built in via TimerOutputs: `MicroMagnetic.reset_timer()` before a
  run, then `println(MicroMagnetic.timer)` afterwards. `run_step`,
  `run_until`, and per-term field evaluation already carry `@timeit` sections.
- **GPU timings are not wall-clock truth**: explicit synchronisation was removed
  from the hot path, so per-section times on GPU backends are lower bounds with
  unknown overlap. Use profiling tools (`nsys`, `ncu` on CUDA) for real numbers.
- Discipline for any optimization: record a baseline (the mumax3 standard
  problems, a 3D FFT-heavy case, a 2D skyrmion case, NEB) → change one thing →
  re-measure.

## Kernel rules (invariants, not suggestions)

1. **No FP64 literals in Float32 kernels.** `1.0 / x::Float32` promotes the
   whole expression to Float64; on GeForce FP64 throughput is 1/32 of FP32, so a
   handful of promoted tail operations cost as much as the entire neighbor
   loop. Write `T(1.0)` / `one(T)` / precomputed `inv(T(v))`. The 25 residual
   sites were fixed in `fb6ead9c`; keep new kernels clean.
2. **`@Const` only on never-Fill arguments** (see
   [Parameter System](params.md#kernel-contracts)).
3. **Derive the kernel backend from data, never from the global**
   (`KernelAbstractions.get_backend(sim.spin)`), and never from host-side
   temporaries — a CPU kernel compiled against device arrays fails with scalar
   indexing.
4. Uniform parameters arrive as `Fill` and are **isbits pass-by-value**: do not
   "unpack" them into arrays for convenience.

## Known structure of the cost

- **Per-term launch + accumulate**: `effective_field` launches one kernel per
  interaction and then a `vector_add` per term. Many-terms simulations pay
  launch overhead; merged-reduction and exch+DMI loop fusion (they share the
  neighbor table; DMI's 4 neighbors are a subset of exchange's 6) are the next
  planned wins, designed as launcher-internal fusions so each term keeps its
  separate energy output.
- **Kernel objects are re-instantiated each step** (closure over backend +
  groupsize); KA caches compilations but instantiation is not free. Caching per
  `(backend, groupsize)` is a candidate once profiling shows it matters.
- **Demag FFT**: plans are cached and reused. On CUDA, the implementation
  assumes cuFFT work on the default stream interacts correctly with KA kernels;
  other GPU backends must re-validate this assumption before trusting timings.
- **Stencil partition layer** (`micro/stencil.jl`): class/pair tables remove
  divisions and parameter gathers from the exchange/DMI inner loops. Reference
  numbers from the implementation branch: 256² CPU/Float64 inline 569 μs →
  partitioned 499 μs (vs 440 μs for the old pre-Fill scalar kernel); on GPU,
  launch overhead dominates at 256² (~20 μs) and 1024² dense cases were ~1.5×
  faster. Only static DMI partitions; the inline kernel is the semantic
  fallback everywhere else.
- **Stopping criterion**: `relax` computes the full `dm/dt` maximum every step
  just for convergence checking; down-sampling is a possible future change (it
  would alter stop timing, so it needs a documented default change).
- **FSAL**: DOPRI54 currently discards its 7th stage, which mathematically
  equals the next step's `k1` — enabling FSAL is a −14% step-time, identity
  transformation (validate against the Butcher table first).

## Saving and I/O

`save_ovf` / `save_vtk` must copy device arrays to the host (`Array(...)`) —
unavoidable, but avoid doing it twice for the same data on hot paths (e.g. when
one step both stores a snapshot and contributes to an energy column).
