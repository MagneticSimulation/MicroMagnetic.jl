# JuMag.jl

_A Julia package for classical spin dynamics and micromagnetic simulations with GPU support._

[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://ww1g11.github.io/JuMag.jl/latest/)
[![Build Status](https://travis-ci.org/ww1g11/JuMag.jl.svg?branch=master)](https://travis-ci.org/ww1g11/JuMag.jl)
[![pipeline status](https://gitlab.com/JuliaGPU/JuMag.jl/badges/master/pipeline.svg)](https://gitlab.com/JuliaGPU/JuMag.jl/commits/master)
[![codecov](https://codecov.io/gl/ww1g11/JuMag.jl/branch/master/graph/badge.svg)](https://codecov.io/gl/ww1g11/JuMag.jl)



### Features

- Support classical spin dynamics and micromagnetic simulations.
- Easily switch between GPU and CPU.
- Easily switch between single and double using `JuMag.cuda_using_double(false)`
- ...

## Installation

Requirements:

- Julia 1.0 (or above) (<http://julialang.org/downloads/>)
- FFTW
- WriteVTK

To enable the GPU support, packages related to CUDA are needed as well:

- CuArrays (optional for GPU support)
- CUDAnative (optional for GPU support)

In [Julia](http://julialang.org), packages can be easily installed using

```
using Pkg;
Pkg.add("FFTW")
```

To install [JuMag.jl](https://github.com/ww1g11/JuMag.jl), download it and append the parent path of JuMag.jl to `JULIA_LOAD_PATH`

```
export JULIA_LOAD_PATH=/parent/path/of/JuMag:$JULIA_LOAD_PATH
```

Now we will see similar messages if we type `using JuMag`

```
julia> using JuMag
[ Info: Precompiling JuMag [8b6b6816-cea2-582c-a99f-83810c20db0f]
┌ Warning: CUDA is not available!
└ @ JuMag ~/Softwares/JuMag.jl/src/JuMag.jl:41
```

Please note in this case `JULIA_LOAD_PATH` is set to `~/Softwares`.
