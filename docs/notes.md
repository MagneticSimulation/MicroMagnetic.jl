Reducing the startup time
-------------------------------

Julia is a dynamically-typed langangue, so the input script will be compiled when we start a simulation. However, the typical startup time in our case ranges from 3s to 30s depends on the complexity of the problem. It is painful especially if we run the simulation using GPU. Luckily, we can compile our package using [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl):

```
using PackageCompiler
compile_package("JuMag")
```

After finishing the compilation, a `dyn.so` file will be generated. If we start julia using `julia -J /path/to/dyn.so` the stratup time will be ignorable.

Note: If you got an error similar to that shown at https://github.com/JuliaLang/PackageCompiler.jl/issues/184, using `dev PackageCompiler` can solve the issue.


Tutorial
-------------------------------
To start a micromagnetic simulation, we first create a FDMesh

```julia
mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100)
```
After that, we create a simulation

```julia
sim = Sim(mesh, name="vortex")
```
and set the damping to 0.5 and switch off the precession term in LLG equation:
```julia
sim.driver.alpha = 0.5
sim.driver.precession = false
```

The geometry of the system can be defined by

```julia
set_Ms(sim, circular_Ms)
```
where `circular_Ms` could be a scalar or a function. The function should take six parameters `(i,j,k,dx,dy,dz)`,
for instance

```julia
function circular_Ms(i,j,k,dx,dy,dz)
	if (i-50.5)^2 + (j-50.5)^2 <= 50^2
		return 8.6e5
	end
	return 0.0
end
```

We add the exchange interaction and the demagnetization field to the system.
```julia
add_exch(sim, 1.3e-11)
add_demag(sim)
```

We need to initialise the system which can be done by defining a function

```julia
function init_fun(i,j,k,dx,dy,dz)
  x = i-50.5
  y = j-50.5
  r = (x^2+y^2)^0.5
  if r<5
    return (0,0,1)
  end
  return (x/r, -y/r, 0)
end
```
and using
```julia
init_m0(sim, init_fun)
```

To trigger the simulation we relax the system

```julia
relax(sim, maxsteps=1000)
```
