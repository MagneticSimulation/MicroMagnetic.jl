# JuMag.jl

_Some Julia codes for classical spin dynamics and finite difference micromagnetics._

[![Build Status](https://travis-ci.org/ww1g11/JuMag.jl.svg?branch=master)](https://travis-ci.org/ww1g11/JuMag.jl) [![Coverage Status](https://coveralls.io/repos/github/ww1g11/JuMag.jl/badge.svg?branch=master)](https://coveralls.io/github/ww1g11/JuMag.jl?branch=master)

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
