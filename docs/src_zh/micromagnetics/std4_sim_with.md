```@meta
ShareDefaultModule = true
```

# 标准问题 4（sim_with）

链接：<https://www.ctcms.nist.gov/~rdm/std4/spec4.html/>


本示例演示如何使用 **MicroMagnetic.jl** 模拟标准问题 4（薄膜几何）。
我们在 GPU 上运行模拟、处理结果并可视化磁化动力学。

导入 MicroMagnetic.jl

````@example
using MicroMagnetic
````

启用 GPU 加速

````@example
@using_gpu()
````

定义体系几何：膜厚 t = 3 nm，长 L = 500 nm，宽 d = 125 nm。
汇总标准问题 4 的所有参数：

````@example
args = (name="std4", task_s=["relax", "dynamics"],           # 任务列表
        mesh=FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9), Ms=8e5, # 饱和磁化强度
        A=1.3e-11,                              # 交换常数
        demag=true,                             # 启用退磁场
        m0=(1, 0.25, 0.1),                      # 初始磁化
        alpha=0.02,                             # Gilbert 阻尼
        steps=100,                              # 动力学步数
        dt=0.01ns,                              # 时间步长
        stopping_dmdt=0.01,                     # 弛豫停止判据
        dynamic_m_every=1,                      # 每步保存磁化
        H_s=[(0, 0, 0), (-24.6mT, 4.3mT, 0)]);
nothing #hide
````

用 `sim_with` 函数运行模拟：

````@example
sim = sim_with(args);
nothing #hide
````

至此标准问题 4 的模拟完成。接下来处理数据，
例如可视化磁化分布，或用模拟结果生成动画。

用 CairoMakie 绘制磁化随时间的演化
````@example
using CairoMakie
fig = plot_ts("std4_llg.txt", ["m_x", "m_y", "m_z"];  xlabel="Time (ns)", ylabel="m", transparency=true)
````

从文件夹 `std4_LLG` 中的模拟 ovf 文件生成动画：
````@example
ovf2movie("std4_LLG"; output="../public/std4.mp4", component='x');
nothing #hide
````
![](../public/std4.mp4)

```@setup
save("../public/std4.png", fig)
```
