# API -- 数据与工具

```@meta
CurrentModule = MicroMagnetic
```

## 数据保存

```@docs
#save_m
save_vtk
save_ovf
read_ovf
mag2ovf
read_table
```


## 工具/可视化

```@docs
voronoi
ovf2vtk
plot_ts
plot_m
ovf2png
ovf2movie
compute_skyrmion_number
compute_guiding_center
compute_magnetic_phase
```


## 服务/GUI

```@docs
start_server
gui
```

## LTEM 工具

```@docs
compute_magnetic_phase
LTEM
MicroMagnetic.magnetic_phase_fft
MicroMagnetic.rotate3d
MicroMagnetic.warp2d
MicroMagnetic.euler_matrix
MicroMagnetic.project3d
MicroMagnetic.vector_field_projection
MicroMagnetic.vector_padding
MicroMagnetic.pad_array
MicroMagnetic.rotation_grid_size
MicroMagnetic.electron_wavelength
MicroMagnetic.interaction_constant
MicroMagnetic.compute_electric_phase
```
