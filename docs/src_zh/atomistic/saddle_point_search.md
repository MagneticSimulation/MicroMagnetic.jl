```@meta
ShareDefaultModule = true
```

# 鞍点搜索（SPS）

系统性的鞍点搜索（SPS）从一个弛豫好的磁化组态出发。它沿低阶 Hessian 模式走到一阶鞍点，
对每个鞍点两侧分别弛豫，并用测地 NEB（GNEB）方法精细化得到的连接。过渡机制从最终的
组态和路径中辨识，而不是预先分配给某个 Hessian 模式。

本例使用 Müller *et al.* 的原子斯格明子模型：开放 ``40 × 40 × 1`` 格点、``J = 1 meV``、
``D/J = 0.45``、``B_z = 0.8 D²/(μ_s J)``。完整可运行脚本见
[`examples/skyrmion_sps.jl`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/examples/skyrmion_sps.jl)。

## 定义模型

```julia
using MicroMagnetic
using CairoMakie

@using_gpu()

const N = 40
const J = 1.0meV
const D = 0.45J
const MU_S = mu_B
const B_Z = 0.8D^2 / (MU_S * J)

function build_sim()
    mesh = CubicMesh(; nx=N, ny=N, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", save_data=false)
    set_mu_s(sim, MU_S)
    add_exch(sim, J)
    add_dmi(sim, D)
    add_zeeman(sim, (0, 0, B_Z))
    return sim
end
```

在开始 SPS 之前，先初始化并弛豫一个斯格明子。示例脚本中包含下文用到的
`initial_skyrmion` 函数。

```julia
sim = build_sim()
init_m0(sim, initial_skyrmion)
relax(sim; max_steps=50_000, stopping_dmdt=1e-5,
      save_data_every=-1, save_m_every=-1)
minimum = Float64.(Array(sim.spin))
```

## 搜索并绘制过渡路径

```julia
result = find_transitions(
    build_sim, minimum;
    n_modes=8,
    directions=(-1, 1),
    exploration_depth=1,
    n_transitions=3,
    images=21,
    output_folder="skyrmion_sps",
)

plot_transition_paths(
    result; output="skyrmion_sps/transition_paths.png")
```

`find_transitions` 会记录失败的候选并继续探索剩余的模式/方向组合。一个被接受的
路径具有两个不同的端点极小、一个一阶鞍点，以及收敛的普通与攀爬像 GNEB 链。能量与力
内部以焦耳为单位；路径图默认以 meV 显示相对能量。

输出根目录包含 `minimum.ovf`、`hessian_modes.txt`、`attempts.txt` 和
`transitions_energy.txt`。每个被接受的 `transition_XX` 文件夹包含鞍点、21 个路径
image 的 OVF 文件、`energy.txt` 和 `distance.txt`；被拒绝的候选不会创建文件夹。与旧版
`NEB` 输出一致，能量表存储以焦耳为单位的绝对能量，距离表存储无量纲测地距离。合并表
额外给出相对路径极小的能量（焦耳与 meV），绘图使用 ``(E-E_min)/meV``，与当前 NEB
文档示例一致。仅在需要 MMF 与 GNEB 迭代历史时才设置 `save_diagnostics=true`。
可在 `@using_gpu()` 之前加载 CUDA、AMDGPU、oneAPI 或 Metal 软件包；否则脚本在 CPU
上运行。

!!! note
    完整的 ``40 × 40`` 搜索属于生产级示例，不会在每次文档构建时执行。其数值结果
    取决于物理模型与容差，但用户通常只需修改模型参数和 `find_transitions` 的参数，
    无需改动求解策略。

## 参考文献

- G. P. Müller, P. F. Bessarab, S. M. Vlasov, F. Lux, N. S. Kiselev,
  S. Blügel, V. M. Uzdin, and H. Jónsson, "Duplication, Collapse, and Escape
  of Magnetic Skyrmions Revealed Using a Systematic Saddle Point Search
  Method," *Physical Review Letters* **121**, 197202 (2018).
  [doi:10.1103/PhysRevLett.121.197202](https://doi.org/10.1103/PhysRevLett.121.197202)

- P. F. Bessarab, V. M. Uzdin, and H. Jónsson, "Method for finding mechanism
  and activation energy of magnetic transitions, applied to skyrmion and
  antivortex annihilation," *Computer Physics Communications* **196**,
  335-347 (2015).
  [doi:10.1016/j.cpc.2015.07.001](https://doi.org/10.1016/j.cpc.2015.07.001)

- G. Henkelman and H. Jónsson, "Improved tangent estimate in the nudged
  elastic band method for finding minimum energy paths and saddle points,"
  *Journal of Chemical Physics* **113**, 9978-9985 (2000).
  [doi:10.1063/1.1323224](https://doi.org/10.1063/1.1323224)

- G. P. Müller *et al.*, "Spirit: Multifunctional framework for atomistic
  spin simulations," *Physical Review B* **99**, 224414 (2019).
  [doi:10.1103/PhysRevB.99.224414](https://doi.org/10.1103/PhysRevB.99.224414)

- H. Schrautzer, M. Sallermann, P. F. Bessarab, and H. Jónsson,
  "Identification of mechanisms of magnetic transitions using an efficient
  method for converging on first-order saddle points," *Physical Review B*
  **112**, 104433 (2025).
  [doi:10.1103/z673-hhnp](https://doi.org/10.1103/z673-hhnp)
