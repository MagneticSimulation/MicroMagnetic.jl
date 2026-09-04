# 数据输入/输出

MicroMagnetic 以 **OVF** 格式(OOMMF 与 mumax3 使用的矢量场格式)写出磁化快照,
以 **VTK** 写出网格与场(供 ParaView 查看),以纯文本**表格**写出标量时间序列。
本页列出可用的读写接口。

## OVF 文件

OVF 是推荐的交换格式:可以与 OOMMF、mumax3 的数据互通,[标准问题](micromagnetics/std4.md)
教程也使用它。

```julia
# 写:当前磁化快照
save_ovf(sim, "m0.ovf")                     # 数据存为 Float64
save_ovf(sim, "m0.ovf"; type=Float32)       # 更小的文件

# 读:既可以直接读进模拟,也可以读进容器
read_ovf(sim, "m0.ovf")                     # 把磁化载入 sim.spin
ovf = read_ovf("m0.ovf")                    # 返回 OVF2 结构体
nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
m = ovf.data                                # 扁平 3N 向量,按 [mx,my,mz,...] 交错
```

`OVF2` 对象可用 `save_ovf(ovf, "copy.ovf")` 写回。要把自己后处理得到的扁平数组
转成 OVF 文件,用 [`mag2ovf`](@ref):

```julia
mag2ovf(m, nx, ny, nz; dx=2.5e-9, dy=2.5e-9, dz=3e-9, fname="processed")
```

`read_ovf` 得到的容器也是分析工具的输入——`compute_skyrmion_number("m.ovf")`
等用法见[分析工具](tools.md)页。

## VTK 文件

VTK 用于 ParaView 风格的网格、磁化与场可视化:

```julia
save_vtk(sim, "config")            # 网格 + m + Ms(开启 demag 电荷时一并保存)
save_vtk(mesh, shape, "shape1")    # 在网格上可视化一个 Shape(见 basics 页)
```

已有的 OVF 快照可以批量转换为 ParaView 可读格式:

```julia
ovf2vtk("skx_LLG")                 # 转换文件夹内所有 OVF
```

## 时间序列表格

每次运行都会向纯文本表 `"<name>_<driver>.txt"` 追加时间、平均磁化与各能量项
(样例见[基础页](basics.md#数据表))。读回为字典:

```julia
data, units = read_table("std4_llg.txt")
data["time"]; data["m_x"]; data["E_total"]
```

运行中追加自定义列用 `SaverItem`——配方(引导中心示例)见
[基础页](basics.md#数据表)。

## 影片与画图

```julia
ovf2movie("std4_LLG"; output="std4.mp4", component='x')  # OVF 文件夹 -> mp4
ovf2png("std4_LLG")                                      # OVF 文件夹 -> png
plot_m(sim; component='z')                               # 实时查看 sim
plot_ts("std4_llg.txt")                                  # 快速画表
```

## 什么时候用什么格式?

| 场景 | 格式 |
| :--- | :--- |
| 与 OOMMF/mumax3 交换快照、喂给分析工具 | OVF |
| ParaView 三维可视化、网格、矢量场 | VTK |
| 在 Julia/Python 里读的时间序列 | 文本表格 |
| 论文/报告用影片 | `ovf2movie` 生成 mp4 |
