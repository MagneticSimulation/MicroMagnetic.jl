# JuMag.jl

_A Julia package for classical spin dynamics and micromagnetic simulations with GPU support._

[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://ww1g11.github.io/JuMag.jl/)
[![Build Status](https://app.travis-ci.com/ww1g11/JuMag.jl.svg?branch=master)](https://app.travis-ci.com/ww1g11/JuMag.jl)
[![Actions Status](https://github.com/ww1g11/JuMag.jl/workflows/CI/badge.svg)](https://github.com/ww1g11/JuMag.jl/actions)
[![codecov](https://codecov.io/github/ww1g11/JuMag.jl/branch/master/graph/badge.svg?token=2t4oGYcWUu)](https://codecov.io/github/ww1g11/JuMag.jl)


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

In [Julia](http://julialang.org), packages can be easily installed using

```
using Pkg;
Pkg.add("CUDA")
```
or

```
julia> ]
(v1.9) pkg> add CUDA
```

To enable GPU support, one has to install one of the following packages:

 | GPU Manufacturer      | Julia Package           |
 | :------------------- | :------------------------ |
 | NVIDIA         | [NVIDIA CUDA](https://github.com/JuliaGPU/CUDA.jl)    |
 | AMD              | [AMD ROCm](https://github.com/JuliaGPU/AMDGPU.jl)    |
 | Intel            | [Intel oneAPI](https://github.com/JuliaGPU/oneAPI.jl) |
 | Apple            | [Apple Metal](https://github.com/JuliaGPU/Metal.jl)  |


To install [JuMag.jl](https://github.com/ww1g11/JuMag.jl), simply using

```
(v1.9) pkg> add https://github.com/ww1g11/JuMag.jl
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

geo = Cylinder(radius=50e-9, height=40e-9) #Create a cylindrical shape with a diameter of 100 nm and a height of 40 nm
mesh = FDMesh(nx=80, ny=80, nz=30, dx=2e-9, dy=2e-9, dz=2e-9) #Create a finite difference mesh to trigger the simulation

sim = create_sim(mesh, shape=geo, Ms=3.87e5, A = 8.78e-12, D = 1.58e-3, demag=true) #create a Sim instance with Fe parameters
init_m0_random(sim) #Initialize a random state

relax(sim, maxsteps=5000, stopping_dmdt=0.1) #Relax the system to obtain a stable magnetization distribution
save_vtk(sim, "m_demag", fields=["demag"]) # Save the magnetization and the stray field into vtk.
```
The magnetization and the stray field around the cylindrical sample are stored in `m_demag.vts`, which can be opened using Paraview. 
