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

## 洛伦兹电镜(LTEM)

LTEM 工具可以从微磁磁化组态出发,模拟洛伦兹透射电镜观测到的物理量:磁相位,
以及可选的 Fresnel 散焦像,并支持**倾斜样品**几何。参考文献:Walton 等,
J. Phys. D **46**, 175005 (2013)(MALTS);Keimpema 等(2006)(平均内电势
相位);倾角/投影方法参考 [MagRecon/maglab](https://github.com/MagRecon/maglab)。

[典型磁结构的 LTEM 成像](ltem.md)一节给出了磁涡旋与 Bloch/Néel skyrmion
的 Fresnel 图像完整算例。

### 磁相位

束流沿 `+z` 时,磁相位只依赖沿束流积分的面内磁化 `T = ∫ M⊥ dz`,以线性填充的
FFT 卷积实现(采样核对离散求和严格成立,并与实空间直接求和逐点交叉验证):

```julia
ovf = read_ovf("skx.ovf")
phi = compute_magnetic_phase(ovf; Ms=8e5)                     # 束流沿 z
phi = compute_magnetic_phase(ovf; Ms=8e5, ty=deg2rad(20))     # 绕 y 倾斜
phi = compute_magnetic_phase(m, deg2rad(20), "y"; Ms=8e5)     # 数组方法
```

`Ms` 单位 A/m(默认 `1/mu_0` 给出单位感应强度的相位);返回 `N × N` 的相位,
单位 rad,样品居中。倾角为绕 `x`/`y`/`z` 轴的欧拉角(`tx`、`ty`、`tz`,弧度):
体数据(连同磁化矢量)通过逐轴双线性切片旋转,束流保持不动。若数据使用别的
束流轴,用 `axis`(`:x`、`:y` 或 `:z`,默认 `:z`)声明,会先转到位再施加倾角。
`N` 可覆盖立方旋转/输出网格尺寸(默认取能容纳旋转后包围盒的 2 的幂)。

!!! note "倾角约定"
    投影得到的面内分量按**实验室系**送入相位核。若按 MagRecon/maglab 的做法把
    投影矢量转回样本系,会混入束流平行分量:均匀面内磁化膜的衬度会放大
    ``1/\cos\alpha`` 倍,垂直磁化膜在任何倾角下衬度恰好为零。实验室系约定
    精确复现洛伦兹偏转 ``\int \mathbf{B}_\perp\,\mathrm{d}z`` ——面内膜
    与倾角无关,垂直膜则 ``\propto \tan\alpha``。

### 完整散焦像

[`LTEM`](@ref) 在投影得到的材料足迹上叠加电相位(平均内电势),再用 Fresnel
散焦传递函数对出射波成像:

```julia
phase, image = LTEM("skx.ovf"; Ms=8e5, V=300, df=200, ty=deg2rad(20))
```

其中 `V` 为加速电压(kV),`df` 为散焦量(µm),`V0` 为平均内电势(V),
`alpha` 为束流发散半角(rad)。返回未解缠的磁相位与归一化强度,均为 `N × N`
且样品居中。`df = 0` 时纯相位对象的强度恒为 1——可作自检。

### 辅助函数

`rotate3d`、`warp2d`、`euler_matrix`、`project3d`、
`vector_field_projection`、`vector_padding`、`electron_wavelength`、
`interaction_constant` 与 `compute_electric_phase` 可单独使用(未导出),
见 `api.md` 的 LTEM 一节。

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
