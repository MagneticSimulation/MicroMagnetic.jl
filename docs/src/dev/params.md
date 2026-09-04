# Parameter System: Fill and make_param

```@meta
CurrentModule = MicroMagnetic
```

Every user-facing parameter — `Ms`, `alpha`, `Ku`, `D`, `temperature`, `pins`,
exchange and DMI constants — accepts a `NumberOrArrayOrFunction`. This page
documents how that union is materialised internally, and the invariants the
implementation relies on.

## Design rule

> **Spatiality is data, not type.** A parameter is either a `Fill` (O(1)
> representation of a uniform array) or a dense array. "Uniform" is one legal
> array representation, not a separate type. No `Spatial*` twin types are
> introduced anywhere.

Before the Fill refactor the type system encoded uniform vs spatial
(`UniformExchange` vs `SpatialExchange`, `LLG` vs `SpatialLLG`, ...), which
duplicated kernels, structs and call sites, and materialised O(N) storage for
every uniform parameter. The single remaining discrimination point is the
setter moment: `set_alpha(sim, 0.1)` stores `Fill(0.1, n)`, while a function or
array is materialised into a dense array. After construction the type is fixed —
no per-step dispatch.

## The Fill type

`src/fill.jl` vendors a ~40-line `Fill{T,N} <: AbstractArray{T,N}` instead of
depending on FillArrays.jl, because kernel-argument behaviour must be controlled
exactly (see the `@Const` rules below):

```julia
struct Fill{T,N} <: AbstractArray{T,N}
    value::T
    size::NTuple{N,Int}
end
```

Properties the implementation relies on:

- **isbits when `T` is isbits** — a `Fill{Float32,1}` kernel argument is passed
  by value, `f[i]` compiles to the same machine code as a scalar, and saves a
  device load;
- O(1) `sum` (`value * prod(size)`), used by `average_m`;
- **immutable** — `copyto!(fill, ...)`, `.=`, `.*=` on a Fill throw; all write
  sites must target dense arrays;
- `Fill{AbstractFloat,1}` is *not* isbits (abstract element type): the symbolic
  eigen mode must not pass such Fills to kernels expecting concrete data.

## make_param: the single materialisation point

All setters route through `make_param(T, v, mesh, n_total; scale=1)`
(`src/params.jl`):

- `v isa Number` (and `T` isbits) → `Fill(T(v), n_total)`;
- `v isa Function` or array → materialised on the **host** and pushed to the
  device with a single `copyto!` (GPU arrays cannot broadcast host arrays —
  never sample into a device array directly);
- eltype is always exactly `T` (the precision frozen at `Sim` construction).

`make_vector_param` is the 3-component analogue for vector parameters (`axis`,
DMI directions, ...), returning three dense component arrays.

## Kernel contracts

1. **`@Const` annotations**: only parameters that can never be a `Fill` keep
   `@Const` (spin, field, neighbor tables, ...). Anything that may be passed as
   a Fill (`alpha`, `mu0_Ms`, `Ku`, `D`, `temperature`, `pins`, ...) must NOT be
   annotated — `@Const` on an isbits struct has unverified behaviour across
   backends. When adding a kernel, follow this split.
2. **No FP64 literals in kernels**: `1.0 / x::Float32` promotes to Float64
   arithmetic, and GeForce-class GPUs run FP64 at 1/32 throughput. Write
   `T(1.0)`, `one(T)`, or precompute `inv(T(v))` host-side (the stencil layer's
   `inv_ms`). This was a measured +38% step-time regression once; the full
   cleanup landed in `fb6ead9c`.
3. **Property assignment semantics**: `sim.driver.alpha = 0.1` still works via a
   `setproperty!` hook on `LLG` — the scalar is converted in place to a
   same-length `Fill`, arrays get their eltype matched, and functions throw
   (spatial sampling needs the mesh; use `set_alpha`).

## Stencil partition layer (exchange / DMI)

Neighbor-loop terms get an extra fast path when parameters are piecewise-constant
per region (`src/micro/stencil.jl`):

- `sim.mat_class[I] = mu0_Ms[I] > 0 ? compressed_region_id : 0` — class 0 means
  vacuum and absorbs all masking (a class-0 row/column in every table is zero);
  region 0 (background) is **not** class 0: non-vacuum cells on the background
  region are a legal material class.
- Per term and axis, an `(R+1)²` pair table `A_eff(class_I, class_J)` is
  precomputed **iff every cell of a class shares the same parameter value**. If
  the check passes, the partition kernel runs with one byte-class gather and no
  division or parameter gather in the inner loop; otherwise the inline kernel
  runs unchanged (semantic fallback).
- `_MAX_CLASSES = 16`: more regions degrade to a 0/1 vacuum mask and the inline
  path (keeps the table in L1).
- Invalidation: `FDMesh.layout_version` is bumped by both `set_region` methods
  (checked lazily per term); `set_Ms` drops the caches eagerly.
- Limitation: partitioning only applies to **static** DMI (`ft === _static_time`);
  time-modulated DMI always takes the inline kernel.

## Symbolic mode

`set_precision(AbstractFloat)` lets parameters carry `FlatTerm` symbolic values
for the eigenmode workflow. In this mode `Fill` of parameters falls back to
dense materialisation where isbits-ness is required, and the eigen code paths
keep their `Vector{AbstractFloat}` buffers. Do not "simplify" these to concrete
element types.
