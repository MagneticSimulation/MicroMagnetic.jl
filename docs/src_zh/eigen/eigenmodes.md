```@meta
ShareDefaultModule = true
```

# 本征模分析

MicroMagnetic.jl 提供了计算磁体系在静态平衡态附近本征模的工具。线性化的
Landau–Lifshitz–Gilbert（LLG）方程给出雅可比矩阵 $B$，它把切向扰动
$\delta\bar{m}$ 映射到其时间导数：

```math
\frac{d\bar{m}}{dt} = B \, \bar{m}
```

其中 $\bar{m}$ 是 2N 维切向矢量（每个自旋格点两个自由度，因为 $|\mathbf{m}|=1$
去掉了径向方向）。$B$ 的本征值 $\sigma = \sigma_r + i\sigma_i$ 的物理含义如下：

- **虚部** $\sigma_i$ → 振荡频率 $f = |\sigma_i| / (2\pi)$
- **实部** $\sigma_r$ → 阻尼率 $\Gamma = -\sigma_r / (2\pi)$
  （稳定平衡态下 $\sigma_r < 0$，因此 $\Gamma > 0$）

**最弱阻尼模**是 $|\sigma_r|$ 最小的模式，即衰减最慢的模式。它们主导长时动力学，
也是实验中最容易观测到的模式。在 Arpack 中用 `which=:LR`（实部最大）来靶向它们——
由于 $\sigma_r$ 为负，"实部最大"意味着最接近零，即阻尼最弱。对于最低频模式，用
`which=:SM`（模最小），需要通过稀疏化做 shift-invert。

## 理论

把 Gilbert 形式的 LLG 方程

```math
\frac{d\mathbf{m}_i}{dt} = -\gamma \, \mathbf{m}_i \times \mathbf{H}_i^{\rm eff}
+ \alpha \, \mathbf{m}_i \times \frac{d\mathbf{m}_i}{dt}
```

在平衡态 $\mathbf{m}_0$（满足
$\mathbf{m}_0 \times \mathbf{H}_0 = 0$）附近线性化，对每个格点 $i$ 得到：

```math
\delta\dot{\mathbf{m}}_i = -\frac{\gamma}{1+\alpha^2} \Big[
  \mathbf{m}_0 \times \delta\mathbf{H}_i
  + \alpha \, \mathbf{m}_0 \times (\mathbf{m}_0 \times \delta\mathbf{H}_i)
  + \delta\mathbf{m}_i \times \mathbf{H}_0
  + \alpha \, \delta\mathbf{m}_i \times (\mathbf{m}_0 \times \mathbf{H}_0)
\Big]
```

其中 $\delta\mathbf{H}_i = \sum_j (\partial \mathbf{H}_i / \partial \mathbf{m}_j)\, \delta\mathbf{m}_j$
是扰动引起的场变化。扰动 $\delta\mathbf{m}$ 位于垂直于
$\mathbf{m}_0$ 的局部切平面内，因此用旋转矩阵 $R_i$ 及其逆 $R_i^{-1}$
投影到每个格点的 2D 标架上。

## 精度设置

本征模分析使用符号扰动类型（`Epsilon`/`FlatTerm`）对局域相互作用做精确线性化
（无有限差分误差）。这要求在**构造模拟之前**设置 `AbstractFloat` 精度：

```@example
using MicroMagnetic
using LinearAlgebra
using CairoMakie

MicroMagnetic.set_backend("cpu")
MicroMagnetic.set_precision(AbstractFloat)
nothing #hide
```

## API

### `build_matrix`

```@docs
build_matrix
```

该函数支持三种输出模式：

| 模式        | 关键字            | 返回类型          | 内存   | 用例                         |
| :---------- | :---------------- | :---------------- | :----- | :--------------------------- |
| 稠密        | *(默认)*          | `Matrix{Float64}` | O(4N²) | 小体系、全谱                 |
| 稀疏        | `sparse=true`     | `SparseMatrixCSC` | O(NZ)  | 稀疏相互作用                 |
| 无矩阵      | `matrixfree=true` | `LLGJacOperator`  | O(N)   | 大体系、Arpack               |

### `LLGJacOperator`

```@docs
LLGJacOperator
```

### `build_demag_matrix`

```@docs
build_demag_matrix
```

## 例 1：立方各向异性 FMR

我们计算带立方各向异性的单胞颗粒的铁磁共振（FMR）频率，并与解析解比较。
设置遵循 [`test/eigen/test_cubic.jl`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/test/eigen/test_cubic.jl)。

对立方各向异性 $K_c$ 和饱和磁化强度 $M_s$，沿 \[100] 方向的有效各向异性场给出 FMR 频率：

```math
f_{[100]} = \frac{\gamma}{2\pi} (H + K), \quad K = \frac{4 K_c}{\mu_0 M_s}
```

```@example
function setup_cubic(; m0=(0,0,1), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=2e-9)
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5)
    add_cubic_anis(sim, 2e4)
    add_zeeman(sim, H)
    init_m0(sim, m0)
    return sim
end

function fmr_frequency_100(H0)
    sim = setup_cubic(m0=(1,0,0), H=(H0, 0, 0))
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
    return imag(eigvals(B)[2]) / 1e9 / (2*pi)
end

function analytical_100(H; Kc=2e4, Ms=8e5)
    gamma = 2.21e5
    K = 4*Kc/(mu_0*Ms)
    return gamma*(H .+ K) / 1e9 / (2*pi)
end

H = 1e5
f_num = fmr_frequency_100(H)
f_ana = analytical_100(H)

println("Numerical:  $(round(f_num, digits=4)) GHz")
println("Analytical: $(round(f_ana, digits=4)) GHz")
println("Relative error: $(round(abs(f_num - f_ana)/f_ana, sigdigits=4))")
```

### 频率随场的扫描

扫描外场即可复现
[`figures/eigen/fig2`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/figures/eigen/fig2/fmr.jl)
中所示的 FMR 色散关系：

```@example
Hs = range(0.8e5, 2e5; length=20)

freqs_100 = [fmr_frequency_100(H) for H in Hs]

function fmr_frequency_110(H0)
    sim = setup_cubic(m0=(1,1,0), H=(H0/sqrt(2), H0/sqrt(2), 0))
    B = build_matrix(sim, gamma=2.21e5, alpha=0.0)
    return imag(eigvals(B)[2]) / 1e9 / (2*pi)
end

freqs_110 = [fmr_frequency_110(H) for H in Hs]

function analytical_110(H; Kc=2e4, Ms=8e5)
    gamma = 2.21e5
    K = Kc/(mu_0*Ms)
    return gamma*sqrt.(H.^2 .- 2K.*H .- 8K^2) / 1e9 / (2*pi)
end

fig = Figure(; size=(500, 360), fontsize=18)
ax = Axis(fig[1, 1]; xlabel="μ₀H (T)", ylabel="Frequency (GHz)")

scatter!(ax, Hs.*mu_0, freqs_100; marker=:rect, markersize=10,
         strokewidth=1, color=:white, label="[100]")
scatter!(ax, Hs.*mu_0, freqs_110; markersize=10,
         strokewidth=1, color=:white, label="[110]")

lines!(ax, Hs.*mu_0, analytical_100.(Hs); color=:slateblue1, label="Analytical")
lines!(ax, Hs.*mu_0, analytical_110.(Hs); color=:slateblue1)

axislegend(ax; position=(0.95, 0.05), labelsize=14)
save("eigen_fig2.png", fig)
nothing #hide
```

![](eigen_fig2.png)

## 例 2：基于无矩阵算子与 Arpack 的大规模 FMR

对于稠密 $2N \times 2N$ 雅可比矩阵放不进内存的大网格，传 `matrixfree=true`
得到 `LLGJacOperator`，再用 `Arpack.eigs` 求解。该算子实现了 Arpack 所需的完整
`AbstractMatrix` 接口（`size`、`eltype`、3 参与 5 参 `mul!`、`issymmetric`、
`copy`）。

该工作流复现标准问题 4 的 FMR 模式图样，见
[`figures/eigen/fig4`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/figures/eigen/fig4/fmr_std.jl)：

```julia
using MicroMagnetic
using SparseArrays
using Arpack
using CairoMakie
using Printf

# 第 1 步：用 Float64 精度弛豫到平衡态
MicroMagnetic.set_precision(Float64)

function setup_std4(; m0=(0,0,1))
    mesh = FDMesh(nx=24, ny=24, nz=2, dx=5e-9, dy=5e-9, dz=5e-9)
    v = (1, 0.7, 0)
    H = 80e3 .* v ./ sqrt(sum(v.^2))
    sim = Sim(mesh; driver="SD")
    set_Ms(sim, 8e5)
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    add_zeeman(sim, H)
    init_m0(sim, m0)
    return sim
end

sim = setup_std4()
relax(sim, stopping_dmdt=0.0001)
m0 = Array(sim.spin)

# 第 2 步：以 AbstractFloat 精度重建，用于本征模分析
MicroMagnetic.set_precision(AbstractFloat)
sim = setup_std4(m0=m0)

# 第 3 步：构建无矩阵算子（O(N) 存储）
op = build_matrix(sim, gamma=2.211e5, matrixfree=true)

# 第 4 步：用 Arpack 求本征值
#   which=:SM  → 最低频模式（需要经稀疏化做 shift-invert）
#   which=:LR  → 最弱阻尼模式（衰减最慢、最易观测）；
#               对无矩阵算子直接调用 eigs(op)。
function compute_eigenvalues(op, sim; nev=10, which=:SM)
    if which == :SM
        A = SparseArrays.sparse(op)
        vals, vecs = Arpack.eigs(A; nev=2*nev, which=:SM, tol=1e-6,
                                 maxiter=10*sim.n_total, explicittransform=:none)
    else
        vals, vecs = Arpack.eigs(op; nev=2*nev, which=which, tol=1e-6,
                                 maxiter=10*sim.n_total, explicittransform=:none)
    end
    idx = findall(x -> imag(x) > 0, vals)
    return vals[idx], vecs[:, idx]
end

freqs, vecs = compute_eigenvalues(op, sim, nev=10, which=:SM)

# 第 5 步：从切向本征矢量重构 3D 磁化扰动
function reconstruct_modes(vecs, sim)
    N = sim.n_total
    mesh = sim.mesh
    nkeep = size(vecs, 2)
    evecs = reshape(vecs, (2, N, nkeep))
    m3d = zeros(ComplexF64, 3, N, nkeep)
    m0 = reshape(sim.spin, 3, N)
    for i in 1:N
        R = MicroMagnetic.rotation_matrix(m0[1,i], m0[2,i], m0[3,i])
        for j in 1:nkeep
            m3d[:, i, j] = R * evecs[:, i, j]
        end
    end
    return reshape(m3d, 3, mesh.nx, mesh.ny, mesh.nz, nkeep)
end

m = reconstruct_modes(vecs, sim)
freqs_GHz = imag.(freqs) ./ (2pi * 1e9)

# 绘制前几个模式的 |m_y|
fig = Figure(size=(600, 290), fontsize=14)
for i in 1:min(5, length(freqs_GHz))
    vv = abs.(m[2, :, :, 1, i])
    title = @sprintf("%.3f GHz", freqs_GHz[i])
    ax = Axis(fig[1, i]; title=title, aspect=1)
    hidedecorations!(ax)
    heatmap!(ax, vv; colormap=:coolwarm, interpolate=true)
end
save("eigen_fig4.png", fig)
```

!!! note "无矩阵算子的内部实现"
`LLGJacOperator` 以 $O(N)$ 内存应用雅可比，步骤如下：

```
1. **展开**：用预计算的旋转矩阵 `R_i` 把 2N 维切向矢量 `x` 展开为
   3N 维笛卡尔 `dm`。
2. **局域场变化** `δH_loc`：对 `m0 + ε·dm` 做一次带 `Epsilon` 标签的符号
   pass，提取有效场的 ε 系数，一次扫描得到 `(∂H_loc/∂m)·dm`。
3. **长程场变化** `δH_demag`：每个退磁相互作用调用一次 `effective_field`
   （FFT 加速）。
4. **LLG 叉乘**：`m0×δH`、`dm×H0` 及阻尼项 `α·m0×(m0×δH)` 等，
   使用缓存的基线场 `H0` 和 `m0×H0`。
5. **投影**：用 `R_i^{-1}` 把 3D 结果投影回 2N 维切向标架。

与稠密构造器（调用 `mul!` 2N 次以实体化所有列）相比，无矩阵的 `mul!`
每次应用只需 **1 次局域符号 pass + 1 次 FFT 退磁调用**。
```

!!! warning "之后要恢复精度"
    `set_precision(AbstractFloat)` 修改的是全局设置。由于 `AbstractFloat`
    不是具体类型，此后创建的模拟会失败（GPU 上尤其如此）。本征模分析结束后，
    请恢复默认值：

```@example
MicroMagnetic.set_precision(Float64)
nothing #hide
```

## 参考文献

- X. Fan, S. Zhang, W. Wang, L. Kong, Y. Guo, Y. Liu, and H. Du, "Automatic
  eigenvalue method in micromagnetic and atomistic simulations,"
  *Journal of Magnetism and Magnetic Materials* **622**, 172989 (2025).
  [doi:10.1016/j.jmmm.2025.172989](https://doi.org/10.1016/j.jmmm.2025.172989)
