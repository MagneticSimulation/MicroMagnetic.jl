```@meta
ShareDefaultModule = true
```

# 纳米条

本例研究尺寸为 60nm × 10nm × 5 nm 的纳米条。
我们用它演示如何用 MicroMagnetic.jl 获得磁化分布。
先导入 MicroMagnetic 和用于绘图的 CairoMakie。

````@example
using MicroMagnetic
using CairoMakie
````

创建 FDMesh

````@example
mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);
nothing #hide
````

用 `Sim` 函数创建带 `SD` 驱动的 Sim，并设置饱和磁化强度 Ms

````@example
sim = Sim(mesh; driver="SD")
set_Ms(sim, 8e5)   #设置饱和磁化强度 Ms=8e5 A/m
````

考虑两种能量（交换与退磁），交换常数 A = 1e-12 J/m。

````@example
add_exch(sim, 1e-12);   #添加交换相互作用
add_demag(sim);         #添加退磁场
nothing #hide
````

把磁化初始化到 (1,1,0) 方向，

````@example
init_m0(sim, (1, 1, 0));  #初始化磁化
nothing #hide
````

可以用 `plot_m` 函数绘制磁化

````@example
plot_m(sim)
````

弛豫体系得到磁化分布。停止判据是 stopping_dmdt，
其取值一般在 [0.01, 1] 范围内。

````@example
relax(sim; max_steps=2000, stopping_dmdt=0.01)
````

再次绘制磁化

````@example
fig = plot_m(sim)
````

把图像保存为 png。

```julia
save("bar.png", fig)
```

保存磁化组态供后续后处理，可用 Paraview（https://www.paraview.org/）可视化

```julia
save_vtk(sim, "bar"; fields=["exch", "demag"])
```

## 使用 sim_with 函数

可以用 sim_with 简化模拟的搭建。
把所有参数放在一起：

````@example
args = (
    task = "Relax",
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2),
    Ms = 8e5,
    A = 1e-12,
    demag = true,
    m0 = (1, 1, 0),
    stopping_dmdt = 0.01
);
nothing #hide
````

然后调用 sim_with 函数

````@example
sim = sim_with(args);
nothing #hide
````

绘制磁化

````@example
plot_m(sim)
````
