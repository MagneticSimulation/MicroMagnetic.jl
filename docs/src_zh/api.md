# API -- 模拟与场

> 本页的分组与导航为中文；各函数的详细文档由源码 docstring 自动生成，**内容为英文**。
> 常用函数的一句话中文语义速查页后续补充。

```@meta
CurrentModule = MicroMagnetic
```

## 网格 -- 微磁模型
```@docs
FDMesh()
```

## 网格 -- 原子模型
```@docs
CubicMesh()
TriangularMesh()
CylindricalTubeMesh()
```

## 形状
```@docs
Plane
Sphere
Cylinder
Box
Torus
```

## 磁化初始化
```@docs
init_m0_random
vortex
skyrmion
bubble2d
skyrmion_lattice
hopfion
```

## 接口
```@docs
set_region
region_map
sim_with
set_backend
set_precision
set_verbose
Sim
NEB
set_driver
set_alpha(sim::AbstractSim, alpha::ArrayOrFunction)
set_pinning(sim::AbstractSim, ids::ArrayOrFunction)
init_m0
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
relax
run_sim
hysteresis
```

## 接口 -- 微磁模型
```@docs
set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)
set_Ms(sim::AbstractSim, shape::Shape, Ms::Number)
add_exch(sim::MicroSim, A::NumberOrTupleOrArrayOrFunction; name="exch")
add_twin_mono_anis
add_magnetoelastic
add_dmi(sim::MicroSim, D::NumberOrTupleOrArrayOrFunction; name="dmi", type="bulk")
add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0, fft=true)
add_exch_int(sim::MicroSimFE, J::Number; k1=1, k2=2, name="exch_int")
add_exch_int(sim::MicroSim, J::NumberOrArrayOrFunction; k1=1, k2=-1, name="exch_int")
add_dmi_int(sim::MicroSim, D::Tuple{Real,Real,Real}; k1=1, k2=-1, name="dmi_int")
rm_demag_charges(sim::MicroSim, Ms; x::Tuple=(0, 0), y::Tuple=(0, 0), z::Tuple=(0, 0), name = "demagc")
voronoi(mesh; min_dist = 20, seed=123456)
```

## 接口 -- 原子模型
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

## 蒙特卡罗
```@docs
run_mc
```

## 其他

```@docs
MicroMagnetic.MicroSim
MicroMagnetic.AtomisticSim
NumberOrArrayOrFunction
NumberOrTupleOrArrayOrFunction
TupleOrArrayOrFunction
```
