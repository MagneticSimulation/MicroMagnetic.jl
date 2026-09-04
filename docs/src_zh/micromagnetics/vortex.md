```@meta
ShareDefaultModule = true
```

# 磁涡旋

导入 MicroMagnetic，并在模拟中使用双精度浮点

````@example
using MicroMagnetic
using CairoMakie
````

本例研究的体系是圆形纳米盘。由于使用有限差分方法，
我们创建尺寸为 200nm × 200nm × 20 nm 的 FDMesh。

````@example
mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=5e-9, nx=100, ny=100, nz=4);
nothing #hide
````

创建直径 100 nm 的圆柱形状

````@example
geo = Cylinder(; radius=100e-9);
nothing #hide
````

开始模拟前需要给定初始态。
这里用函数给出初始态，函数接受六个参数 `(i,j,k,dx,dy,dz)`。

````@example
function init_fun(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end
````

我们定义

````@example
args = (
    task = "Relax",
    mesh = mesh,
    shape = geo,
    Ms = 8e5,
    A = 1.3e-11,
    demag = true,
    m0 = init_fun,
    stopping_dmdt = 0.01
);

sim = sim_with(args);
nothing #hide
````

用 plot_m 函数绘制磁化分布

````@example
plot_m(sim; figsize=(400, 400), arrows=(30, 30))
````

同时保存 vtk

```julia
save_vtk(sim, "vortex")
```

