# References

### Meshes

```@docs
FDMesh
FEMesh
FDMeshGPU
CubicMeshGPU
TriangularMeshGPU
```


### DataTypes

```@docs
JuMag.NumberOrArrayOrFunction
JuMag.NumberOrArray
JuMag.TupleOrArrayOrFunction
```


### Interfaces

```@docs
create_sim
Sim
set_Ms
set_Ms_cylindrical
set_driver
create_box
create_cylinder
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
add_exch_vector
add_exch_kagome
add_anis_kagome
add_magnetoelectric_laser
add_exch_anis
update_zeeman
update_anis
relax
compute_guiding_center
```


### DataSaving

```@docs
save_m
save_vtk
save_ovf
read_ovf
```



### Tools

```@docs
ovf2vtk
OVF2LTEM
OVF2MFM
OVF2XRAY
plot_ovf_slice
plot_ovf_projection
```

