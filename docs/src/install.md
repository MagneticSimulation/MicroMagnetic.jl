## Installation

Install MicroMagnetic is straightforward as long as Julia (<http://julialang.org/downloads/>) is installed, and it is equally easy in Windows, Linux and Mac.  

In [Julia](http://julialang.org), packages can be easily installed with the Julia package manager.
From the Julia REPL, type ] to enter the Pkg REPL mode and run:

```julia
pkg> add MicroMagnetic
```

Or, equivalently:

```julia
julia> using Pkg;
julia> Pkg.add("MicroMagnetic")
```

To install the latest development version:
```julia
pkg> add MicroMagnetic#master
```

To enable GPU support, one has to install one of the following packages:

!!! note "GPU Support"
    | GPU Manufacturer      | Julia Package                                      |
    | :------------------:  | :-----------------------------------------------:  |
    | NVIDIA                | [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)     |
    | AMD                   | [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl) |
    | Intel                 | [oneAPI.jl](https://github.com/JuliaGPU/oneAPI.jl) |
    | Apple                 | [Metal.jl](https://github.com/JuliaGPU/Metal.jl)   |

For example, we can install `CUDA` for NVIDIA GPUs:

```julia
pkg> add CUDA
```

Now we will see similar messages if we type `using MicroMagnetic`

```
julia> using MicroMagnetic
julia> using CUDA
Precompiling CUDAExt
  1 dependency successfully precompiled in 8 seconds. 383 already precompiled.
[ Info: Switch the backend to CUDA.CUDAKernels.CUDABackend(false, false)
```


## Running MicroMagnetic.jl from Python

Thanks to [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl), running MicroMagnetic.jl from Python is seamless.

Below is an example setup for **Standard Problem #4**:

```python
from juliacall import Main as jl
jl.seval("using MicroMagnetic")

# Define simulation parameters
args = {
    "name": "std4",
    "task_s": ["relax", "dynamics"],                   # List of tasks to perform
    "mesh": jl.FDMesh(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9),  # Julia FDMesh object
    "Ms": 8e5,                                         # Saturation magnetization (A/m)
    "A": 1.3e-11,                                      # Exchange stiffness constant (J/m)
    "demag": True,                                     # Enable demagnetization
    "m0": (1, 0.25, 0.1),                              # Initial magnetization vector
    "alpha": 0.02,                                     # Gilbert damping coefficient
    "steps": 100,                                      # Number of dynamic simulation steps
    "dt": 0.01 * jl.ns,                                # Time step size (0.01 ns)
    "stopping_dmdt": 0.01,                             # Stopping criterion for relaxation
    "dynamic_m_interval": 1,                           # Save magnetization at each step
    "H_s": [(0, 0, 0), (-24.6 * jl.mT, 4.3 * jl.mT, 0)]  # Sequence of applied magnetic fields
}

# Run the simulation
sim = jl.sim_with(**args)
```