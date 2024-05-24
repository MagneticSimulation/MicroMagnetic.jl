# API

```@meta
CurrentModule = MicroMagnetic
```

## Meshes
```@docs
FDMesh
CubicMesh
TriangularMesh
CylindricalTubeMesh
```

## Shapes
```@docs
Cylinder
Box
Torus
```

## Interfaces

```@docs
set_backend
create_sim
run_sim
Sim
set_Ms
set_driver
set_mu_s
init_m0
init_m0_random
init_m0_skyrmion
add_exch
add_anis
add_cubic_anis
add_dmi
add_demag
add_zeeman
add_exch_int
add_dmi_int
#add_exch_kagome
#add_anis_kagome
#add_magnetoelectric_laser
#add_exch_anis
#add_exch_rkky
#
update_zeeman
update_anis
relax
#compute_guiding_center
```

## DataSaving

```@docs
#save_m
save_vtk
save_ovf
read_ovf
```


## Tools

```@docs
ovf2vtk
#plot_m
#jdl2png
#jdl2movie
```

## DataTypes

```@docs
MicroMagnetic.NumberOrArrayOrFunction
MicroMagnetic.NumberOrTupleOrArrayOrFunction
MicroMagnetic.NumberOrArray
MicroMagnetic.TupleOrArrayOrFunction
```

## Others

```@docs
MicroMagnetic.MicroSim
MicroMagnetic.AtomisticSim
```