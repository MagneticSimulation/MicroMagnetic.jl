# API for users

```@meta
CurrentModule = MicroMagnetic
```

## Meshes -- MicroMagnetic 
```@docs
FDMesh()
```

## Meshes -- Atomistic 
```@docs
CubicMesh()
TriangularMesh()
CylindricalTubeMesh()
```

## Shapes
```@docs
Cylinder
Box
Torus
```

## Interfaces 
```@docs
sim_with
set_backend
set_precision
set_verbose
Sim
NEB
set_driver
set_alpha(sim::AbstractSim, alpha::ArrayOrFunction)
init_m0
init_m0_random
init_m0_skyrmion
add_zeeman
add_anis
add_cubic_anis
add_hex_anis
update_zeeman
update_anis
add_stt
add_sot
add_torque
add_sahe_torque
add_thermal_noise
create_sim
relax
run_sim
#compute_guiding_center
```

## Interfaces -- MicroMagnetic 
```@docs
set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)
set_Ms(sim::AbstractSim, shape::Shape, Ms::Number)
add_exch(sim::MicroSim, A::NumberOrTupleOrArrayOrFunction; name="exch")
add_dmi(sim::MicroSim, D::NumberOrTupleOrArrayOrFunction; name="dmi", type="bulk")
add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0, fft=true)
add_exch_int(sim::MicroSim, J::Float64; k1=1, k2=-1, name="exch_int")
add_dmi_int(sim::MicroSim, D::Tuple{Real,Real,Real}; k1=1, k2=-1, name="dmi_int")
voronoi(mesh; min_dist = 20, seed=123456)
```

## Interfaces -- Atomistic 
```@docs
set_mu_s(sim::AtomisticSim, init::NumberOrArrayOrFunction)
add_exch(sim::AtomisticSim, J1::NumberOrArray; name="exch", J2=0, J3=0, J4=0)
add_exch(sim::AtomisticSim, Jfun::Function; name="exch")
add_exch_bq(sim::AtomisticSim, K::NumberOrArray; name="exch_bq")
add_dmi(sim::AtomisticSim, D::Real; name="dmi", type="bulk")
add_dmi(sim::AtomisticSim, Dij::Array{<:Real, 2}; name="dmi")
add_dmi(sim::AtomisticSim, Dfun::Function; name="dmi")
add_demag(sim::AtomisticSim; name="demag", Nx=0, Ny=0, Nz=0)
add_exch_kagome(sim::AtomisticSim, Jxy::Number, Jz::Number; name="exch")
add_anis_kagome(sim::AtomisticSim, Ku::Float64; ax1=(-0.5, -sqrt(3) / 2, 0), ax2=(1, 0, 0), ax3=(-0.5, sqrt(3) / 2, 0), name="anis")
add_anis_tube(sim::AtomisticSim, Ku::Float64; name="anis")
```

## DataSaving

```@docs
#save_m
save_vtk
save_ovf
read_ovf
```


## Tools/Visualization

```@docs
voronoi
ovf2vtk
plot_ts
plot_m
ovf2png
ovf2movie
```


## Others

```@docs
MicroMagnetic.MicroSim
MicroMagnetic.AtomisticSim
```