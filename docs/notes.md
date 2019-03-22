Preparing
-------------------------------

We try to use [Julia](http://julialang.org). Installation and IDE stuff can be found at http://julialang.org/downloads/. The suggested IDE is atom, seems quite nice.

Design
-------------------------------
To start a micromagnetic simulation, we first create a FDMesh

```julia
mesh = FDMesh(dx=1.0, dy=1.0, dz=1.0, nx=100, ny=100, nz=1, unit_length=1e-9)
```
and the geometry of the system can be defined by

```julia
set_Ms(mesh, circular_Ms)
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
