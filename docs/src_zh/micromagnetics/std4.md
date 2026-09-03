```@meta
ShareDefaultModule = true
```

# 标准问题 4

本示例演示如何使用 **MicroMagnetic.jl** 模拟标准问题 4。我们先将体系弛豫到稳定的磁化组态，
再施加外磁场研究其影响。

导入所需模块

````@example
using MicroMagnetic
using CairoMakie
````

启用 GPU 加速

````@example
@using_gpu()
````

定义体系几何：膜厚 t = 3 nm，长 L = 500 nm，宽 d = 125 nm。

````@example
mesh = FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);
nothing #hide
````

## 第一步：弛豫体系

模拟的第一步是弛豫体系，得到通常称为 "S" 态的稳定磁化组态。我们把这一过程封装在 `relax_system` 函数中。

````@example
function relax_system(mesh)
   sim = Sim(mesh; driver="SD", name="std4")

   Ms = 8e5
   A = 1.3e-11

   set_Ms(sim, Ms)   # 设置饱和磁化强度
   add_exch(sim, A)  # 添加交换相互作用
   add_demag(sim)    # 添加退磁场

   init_m0(sim, (1, 0.25, 0.1))  # 初始化磁化
   relax(sim; stopping_dmdt=0.01)  # 弛豫体系

   return sim
end
````

`relax_system` 函数以网格为输入，执行以下步骤：

**模拟初始化：**
函数用给定网格初始化模拟（`sim`），并把驱动设为 "SD"（最速下降），该驱动通常用于弛豫过程。模拟命名为 "std4"，与标准问题 4 保持一致。

**材料参数设置：**
分别通过 `set_Ms` 和 `add_exch` 设置饱和磁化强度 `Ms` 和交换常数 `A`，并通过 `add_demag` 加入退磁场效应。

**初始磁化：**
初始磁化矢量 `m0` 通过 `init_m0` 设置，这里取 `(1, 0.25, 0.1)`。

**弛豫：**
用 `relax` 函数弛豫体系，迭代极小化体系能量，直到磁化变化率（`dm/dt`）低于给定阈值（`stopping_dmdt=0.01`）。

弛豫体系，得到稳定的磁化组态

````@example
sim = relax_system(mesh);
nothing #hide
````

## 第二步：施加外磁场

得到稳定的 "S" 态后，下一步是施加外磁场并观察磁化动力学。
用 `plot_m` 函数绘制磁化分布

````@example
plot_m(sim; component='x')
````

从弛豫好的 "S" 态出发施加外磁场

```julia
set_driver(sim; driver="LLG", alpha=0.02, gamma=2.211e5)
add_zeeman(sim, (-24.6mT, 4.3mT, 0))  # 施加外磁场
run_sim(sim; steps=100, dt=1e-11, save_m_every=1)   # 运行 100 步动力学模拟
```

## 第三步：结果可视化

默认情况下，`run_sim` 会生成 `std4_LLG` 文件夹存放磁化数据。可以用它生成动画来观察磁化动力学。
根据模拟结果生成动画

```julia
ovf2movie("std4_LLG"; output="../public/std4.mp4", component='x');
```
![](../public/std4.mp4)
