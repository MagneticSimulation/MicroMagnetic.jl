# 分析工具

MicroMagnetic 自带一组分析例程,可以直接作用于模拟对象、磁化数组或 OVF 文件。
它们位于 `src/tools/`。

## 拓扑荷

[`compute_skyrmion_number`](@ref) 用 Berg–Lüscher 方法计算斯格明子数 *Q*(拓扑荷)。
它接受扁平磁化数组(`sim.spin` 正是这种格式)或直接接受 OVF 文件:

```julia
Q = compute_skyrmion_number(Array(sim.spin), mesh)   # 当前模拟状态
Q = compute_skyrmion_number("snapshot.ovf")          # 磁盘上的任意 OVF 文件
```

分层变体 `compute_skyrmion_number_layers("snapshot.ovf")` 返回每层的荷——对
需要逐层分析斯格明子数的薄膜体系很有用。

## 引导中心

[`compute_guiding_center`](@ref) 返回涡旋或斯格明子的引导中心 `(Rx, Ry)`
(Thiele 坐标)。它接受模拟对象或扁平磁化数组,并可用关键字限制在一个子盒内:

```julia
Rx, Ry = compute_guiding_center(sim)
Rx, Ry = compute_guiding_center(sim; z=3)   # 限制在第 k=3 层
```

要在整个运行过程中追踪引导中心,把它作为 saver 的附加列,而不是事后处理快照:

```julia
item = SaverItem(("Rx", "Ry"), ("<m>", "<m>"), compute_guiding_center)
push!(sim.saver.items, item)
```

`micromagnetics/skyrmion_stt.md` 示例(斯格明子动力学 STT)就是这样记录涡旋核
轨迹的。

## 洛伦兹电镜(LTEM)相位恢复

[`compute_magnetic_phase`](@ref) 对给定磁化组态计算 LTEM 会测得的磁相位,使用
沿指定轴的傅里叶空间投影(Beleggia/Phatak 型方法)。返回 `Nx × Ny` 的相位矩阵,
单位 rad:

```julia
ovf = read_ovf("skx.ovf")
phi = compute_magnetic_phase(ovf; Ms=8e5, axis=:z)   # OVF2 容器
```

另有数组方法 `compute_magnetic_phase(m, theta, axis; N1, N2, Ms, d)`,接受
`(3, nx, ny, nz)` 形状的 4D 磁化数组,可在投影前用欧拉角 `theta` 倾斜样品
(转轴 `"X"` 或 `"Y"`),投影/傅里叶核的填充用 `N1`/`N2` 控制。

!!! note
    LTEM 例程假设样品厚度均匀(取自网格步长),相位约定遵循洛伦兹显微术常用的
    欠焦强度传递。

## Voronoi 镶嵌(多晶样品)

[`voronoi`](@ref) 对网格生成 Voronoi 镶嵌并给每个 cell 分配 region id——这是
多晶模拟的起点,晶界处可以赋予不同的材料参数:

```julia
mesh = FDMesh(; nx=200, ny=200, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9)
voronoi(mesh; min_dist=20, seed=123456)   # 把晶粒 id 填进 mesh.region
plot_voronoi(mesh)                        # 检查镶嵌结果

# 用 region_map 设置逐晶粒参数
set_Ms(sim, region_map(1 => 8.6e5, 2 => 7.9e5))
add_exch(sim, region_map(1 => 1.3e-11, 2 => 1.0e-11))
```

`min_dist`(单位是 cell 数)控制晶粒中心的最小间距,`seed` 保证镶嵌可复现。
这样生成的 region 与 [`set_region`](@ref) 以及[基础页](basics.md)描述的
region 参数函数完全兼容。
