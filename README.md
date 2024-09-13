# MicroMagnetic.jl 

### **Previous name: JuMag.jl**
_A Julia package for classical spin dynamics and micromagnetic simulations with GPU support._


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ww1g11/MicroMagnetic.jl/gh-pages)
[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://ww1g11.github.io/MicroMagnetic.jl/dev/)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ww1g11.github.io/MicroMagnetic.jl/stable/)
[![Actions Status](https://github.com/ww1g11/MicroMagnetic.jl/workflows/CI/badge.svg)](https://github.com/ww1g11/MicroMagnetic.jl/actions)
[![codecov](https://codecov.io/github/ww1g11/MicroMagnetic.jl/branch/master/graph/badge.svg?token=2t4oGYcWUu)](https://codecov.io/github/ww1g11/MicroMagnetic.jl)
[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1674--1056%2Fad766f-blue)](https://doi.org/10.1088/1674-1056/ad766f)

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

## Run MicroMagnetic.jl in the cloud

You don't have to install anything and you can run MicroMagnetic.jl in the cloud via Binder:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ww1g11/MicroMagnetic.jl/gh-pages)

All the julia scripts and Jupyter Notebooks (in tutorials folder) hosted in the repository can be executed and modified.

## Documentation

The documentation is available:

[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://ww1g11.github.io/MicroMagnetic.jl/dev/)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ww1g11.github.io/MicroMagnetic.jl/stable/)

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


## Quick start
Assuming we have a cylindrical FeGe sample with a diameter of 100 nm and a height of 40 nm, we want to know 
its magnetization distribution and the stray field around it. We can use the following script: 

```julia
using MicroMagnetic
@using_gpu() # Import available GPU packages such as CUDA, AMDGPU, oneAPI, or Metal

# Define simulation parameters
args = (
    task = "Relax", # Specify the type of simulation task
    mesh = FDMesh(nx=80, ny=80, nz=30, dx=2e-9, dy=2e-9, dz=2e-9), # Define the mesh grid
    shape = Cylinder(radius=50e-9, height=40e-9), # Define the shape 
    Ms = 3.87e5, # Set the saturation magnetization (A/m)
    A = 8.78e-12, # Set the exchange stiffness constant (J/m)
    D = 1.58e-3, # Set the Dzyaloshinskii-Moriya interaction constant (J/m^2)
    demag = true, # Enable demagnetization effects in the simulation
    m0 = (1,1,1), # Set the initial magnetization direction
    stopping_dmdt = 0.1 # Set the stopping criterion 
);

# Run the simulation with the specified parameters
sim = sim_with(args); 

# Save the magnetization and the stray field into vtk.
save_vtk(sim, "m_demag", fields=["demag"]) 
```
The magnetization and the stray field around the cylindrical sample are stored in `m.vts`, which can be opened using Paraview. 

## Questions and Contributions

If you have any questions about usage, please join the conversation on our [GitHub Discussions page](https://github.com/ww1g11/MicroMagnetic.jl/discussions). 

We greatly appreciate contributions, feature requests, and suggestions! If you encounter any issues or have ideas for improving MicroMagnetic.jl, feel free to open an issue on our [GitHub Issues page](https://github.com/ww1g11/MicroMagnetic.jl/issues).

## Citation

If you use MicroMagnetic.jl in your research,  please cite the following publication:

```
@article{10.1088/1674-1056/ad766f,
  author={Wang, Weiwei and Lyu, Boyao and Kong, Lingyao and Fangohr, Hans and Du, Haifeng},
  title={MicroMagnetic.jl: A Julia package for micromagnetic and atomistic simulations with GPU support},
  journal={Chinese Physics B},
  url={http://iopscience.iop.org/article/10.1088/1674-1056/ad766f},
  year={2024}
}
```

