Preparing
-------------------------------

We try to use [Julia](http://julialang.org). Installation and IDE stuff can be found at http://julialang.org/downloads/. The suggested IDE is atom, seems quite nice.

We probably will make use of [Sundials](https://github.com/jgoldfar/Sundials.jl), which can be installed through

```julia
Pkg.add("Sundials")
```

Design
-------------------------------
It seems that the OOP in Julia is not fully supported, so we will try to avoid
to use them. For example, instead of creating a mesh using the constructors, we
use a function:

```julia
mesh = create_mesh(dx=1.0, nx=10)
```
