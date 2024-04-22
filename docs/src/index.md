# JuMag.jl

_A Julia package for classical spin dynamics and micromagnetic simulations with GPU support._


### Features

- Supports classical spin dynamics and micromagnetic simulations.
- Compatible with CPU and multiple GPU platforms, including NVIDIA, AMD, Intel, and Apple GPUs.
- Supports both double and single precision.
- Supports Monte Carlo simulations for atomistic models.
- Implements the Nudged-Elastic-Band method for energy barrier computations.
- Supports Spin-transfer torques, including Zhang-Li and Slonczewski models.
- Incorporates various energy terms and thermal fluctuations.
- Supports constructive solid geometry.
- Supports periodic boundary conditions.
- Easily extensible to add new features.

## Installation

Install JuMag is straightforward as long as Julia (<http://julialang.org/downloads/>) is installed, and it is equally easy in Windows, Linux and Mac.  

In [Julia](http://julialang.org), packages can be easily installed with the Julia package manager.
From the Julia REPL, type ] to enter the Pkg REPL mode and run:

```julia
pkg> add https://github.com/ww1g11/JuMag.jl
```

Or, equivalently:

```julia
julia> using Pkg;
julia> Pkg.add("https://github.com/ww1g11/JuMag.jl")
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

Now we will see similar messages if we type `using JuMag`

```
julia> using JuMag
julia> using CUDA
Precompiling CUDAExt
  1 dependency successfully precompiled in 8 seconds. 383 already precompiled.
[ Info: Switch the backend to CUDA.CUDAKernels.CUDABackend(false, false)
```

# Quick start
Assuming we have a cylindrical FeG sample with a diameter of 100 nm and a height of 40 nm, we want to know its magnetization distribution and the stray field around it. 
We can use the following script: 

```julia
using JuMag
@using_gpu() #Import available GPU packages such as CUDA, AMDGPU, oneAPI or Metal

# Create a cylindrical shape with a diameter of 100 nm and a height of 40 nm
geo = Cylinder(radius=50e-9, height=40e-9) 

# Create a finite difference mesh to trigger the simulation
mesh = FDMesh(nx=80, ny=80, nz=30, dx=2e-9, dy=2e-9, dz=2e-9) 

# Create a Sim instance with Fe parameters
sim = create_sim(mesh, shape=geo, Ms=3.87e5, A = 8.78e-12, D = 1.58e-3, demag=true) 

init_m0_random(sim) #Initialize a random state

#Relax the system to obtain a stable magnetization distribution
relax(sim, maxsteps=5000, stopping_dmdt=0.1) 

# Save the magnetization and the stray field into vtk.
save_vtk(sim, "m_demag", fields=["demag"]) 
```
The magnetization and the stray field around the cylindrical sample are stored in `m_demag.vts`, which can be opened using Paraview. 

