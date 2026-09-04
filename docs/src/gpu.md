# GPU & Precision

MicroMagnetic.jl runs on the CPU by default and supports NVIDIA, AMD, Intel and Apple
GPUs through [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl).
This page explains the three session-level settings — **backend**, **precision** and
**group size** — and the one rule you need to remember:

!!! warning "Settings are constructor-time defaults"
    `set_backend`, `set_precision` and `set_groupsize` only affect simulations
    created **after** the call. A `Sim` freezes its precision into its type and
    its backend into its data at construction time; later changes to the global
    settings never touch an existing simulation.

## Backends

Install one of the GPU packages listed on the [installation page](install.md) and load it
through the `@using_gpu` macro, which imports whichever backend package is available
in your environment:

```julia
using MicroMagnetic
@using_gpu()   # loads CUDA, AMDGPU, oneAPI or Metal if installed, else stays on CPU
```

The macro switches the default backend automatically. If you need manual control:

```julia
set_backend("cuda")   # allocations after this call go to the GPU
set_backend("cpu")    # back to CPU
```

Running `set_backend` when simulations already exist prints a warning — the running
simulations keep their data where they are.

## Precision

The default working precision is `Float64`. Single precision usually doubles the
throughput on consumer GPUs (FP64 throughput is 1/32 of FP32 on GeForce cards) and
halves the memory footprint:

```julia
set_precision(Float32)   # simulations created after this line use Float32
sim = Sim(mesh; driver="LLG")
```

Allowed values are `Float64`, `Float32` and `AbstractFloat`:

- `Float64` / `Float32` — concrete precision for production runs.
- `AbstractFloat` — a **symbolic mode** used only by the eigenmode workflow; it lets
  parameters carry symbolic terms. Regular simulations fail in this mode, in
  particular on GPUs (device arrays cannot hold an abstract element type).

!!! note "Set precision before creating simulations"
    Because precision is frozen into the simulation type at construction, call
    `set_precision` **before** `Sim(...)` / `CubicMesh(...)` / `FDMesh(...)`.
    Mixing precisions across simulations in one session is not supported.

### The eigenmode exception

The eigenmode analysis requires the symbolic mode, which changes the global
precision for the **whole session** and is not restored automatically. After
finishing an eigenmode calculation, restore the concrete precision before creating
further simulations:

```julia
# ... eigenmode analysis with set_precision(AbstractFloat) ...
set_precision(Float64)   # restore before the next simulation
```

See the [eigenmode analysis](eigen/eigenmodes.md) page for the full workflow.

## Group size

`set_groupsize` controls the thread-group size used when launching GPU kernels
(default 256). It is a pure performance knob — it never affects results:

```julia
set_groupsize(128)   # tune only if profiling suggests it
```

## Checking what is active

After `using MicroMagnetic` the startup message reports the active backend. In a
script you can verify the frozen precision through the simulation itself:

```julia
sim = Sim(mesh; driver="LLG")
eltype(sim.spin)                    # the working precision of this simulation
```

!!! tip "Benchmark before tuning"
    For anything performance related — group size, precision, backend — measure
    with `MicroMagnetic.timer` before and after; the Timings section on the
    [basics page](basics.md) shows how the report looks.
