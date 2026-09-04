# 参与贡献

### 分步指南

#### 1. **Fork 仓库**
Fork 会在你自己的 GitHub 账号下创建一份 MicroMagnetic.jl 仓库的个人副本。

- 打开 **MicroMagnetic.jl** 仓库：[MicroMagnetic.jl](https://github.com/magneticsimulation/MicroMagnetic.jl)。
- 点击右上角的 **Fork** 按钮，在你的 GitHub 账号中创建仓库副本。

#### 2. **克隆你的 Fork**
Fork 之后，把你的 fork 克隆到本地机器上就可以开始工作了。

- 打开终端（Windows 上用 Git Bash）。
- 使用以下命令克隆仓库（把 `yourusername` 换成你的 GitHub 用户名）：
  ```bash
  git clone https://github.com/yourusername/MicroMagnetic.jl
  cd MicroMagnetic.jl
  ```

#### 3. **（可选）把原仓库设为 upstream**
为了让你的 fork 与原仓库保持同步，可以把原始的 `MicroMagnetic.jl` 仓库添加为名为 `upstream` 的 remote。这一步是可选的，但能让你方便地拉取主仓库的更新。

- 运行以下命令进行设置：
  ```bash
  git remote add upstream https://github.com/magneticsimulation/MicroMagnetic.jl
  git fetch upstream
  ```

如果不设置 `upstream`，也照样可以贡献，只是需要手动管理原仓库的更新。

#### 4. **创建新分支**
在修改代码之前，先创建一个新分支，把你的改动与主分支隔离开。建议使用描述性的分支名，例如 `fix-simulation-bug` 或 `add-feature-x`。

- 创建并切换到新分支：
  ```bash
  git checkout -b my-new-feature
  ```

#### 5. **在 Julia 中激活本地开发**
为了让 Julia 使用你本地的包而不是注册表中的版本，需要激活开发模式。假设你把仓库克隆到了 `/path/to/MicroMagnetic.jl`，使用以下命令：

```julia
  using Pkg
  Pkg.develop(path="/path/to/MicroMagnetic.jl")
```

这样就可以使用本地版本的包进行开发了。

#### 6. **修改代码**
现在可以在新分支上修改代码了，例如添加新功能、修复 bug 或更新文档。

- 用你喜欢的编辑器或 IDE 修改 **MicroMagnetic.jl** 仓库中的源码。

#### 7. **测试你的改动**
推送之前，先运行测试确保改动正常。测试通常位于 `test` 文件夹中。

- 运行测试：
  ```bash
  julia --project test/runtests.jl
  ```

确保所有测试通过后再继续。如果有测试失败，先调试并修复代码。

#### 8. **提交改动**
对改动满意后，就可以提交了。请写一条清晰简洁的提交信息说明你做了什么。

- 暂存改动：
  ```bash
  git add .
  ```
- 提交：
  ```bash
  git commit -m "Add feature X to improve simulation accuracy"
  ```

#### 9. **推送到你的 Fork**
本地提交完成后，把分支推送到你在 GitHub 上的 fork 仓库。

- 推送分支：
  ```bash
  git push origin my-new-feature
  ```

#### 10. **发起 Pull Request**
改动推送到 fork 之后，就可以发起 pull request（PR），把改动合并回主仓库。

- 打开你的 fork 仓库页面。
- 会看到发起 pull request 的按钮，点击它。
- 在 PR 中详细描述你做的改动、为什么需要这些改动以及其他相关信息。
- 提交 PR。

#### 11. **（可选）保持 Fork 最新**
随着原 `MicroMagnetic.jl` 仓库的演进，你需要让 fork 与最新改动保持同步。如果按第 3 步设置了 `upstream`，这一步很简单。

- 拉取原仓库的最新更新：
  ```bash
  git fetch upstream
  ```
- 把更新合并到你的本地分支：
  ```bash
  git merge upstream/main
  ```

### 示例：添加六角各向异性

所有改动见[这个提交](https://github.com/MagneticSimulation/MicroMagnetic.jl/commit/6a9d7ec7f01cacfaa7a908dc286f58a5bcae910b)。

#### 1. **找到能量表达式**
六角各向异性的能量密度为：

```math
E = K_1 \sin^2 \theta + K_2 \sin^4 \theta + K_3 \sin^6 \theta \cos 6\phi
```

其中 $\theta$ 是磁化矢量 $\mathbf{m}$ 与 $c$ 轴（$z$ 轴）的夹角。
$\phi$ 是 $\mathbf{m}$ 在六角面上的投影相对于 $x$ 轴的夹角。

#### 2. **计算有效场**
利用恒等式 $\cos 6x = -\sin^6 x + 15 \cos^2 x \sin^4 x - 15 \cos^4 x \sin^2 x + \cos^6 x$，能量密度可以用 $m_x$、$m_y$、$m_z$ 改写为：

```math
E = K_1 (1 - m_z^2) + K_2 (1 - m_z^2)^2 + K_3 \left( m_x^6 - 15m_x^4 m_y^2 + 15m_x^2 m_y^4 - m_y^6 \right)
```

对应的有效场为：

```math
\mathbf{H}_\mathrm{eff} = -\frac{6K_3}{\mu_0 M_s} \left( m_x^5 - 10m_x^3m_y^2 + 5m_xm_y^4 \right) \mathbf{e}_x 
 -\frac{6K_3}{\mu_0 M_s} \left( -6m_x^4m_y + 10m_x^2m_y^3 - m_y^5 \right) \mathbf{e}_y 
+ \frac{2m_z}{\mu_0 M_s} \left[ K_1 + 2K_2(1 - m_z^2) \right] \mathbf{e}_z
```

#### 3. **在 `src/head.jl` 中定义 struct**
```julia
mutable struct HexagonalAnisotropy{T<:AbstractFloat} <: MicroEnergy
    K1::T
    K2::T 
    K3::T
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end
```

#### 4. **在 `src/micro/add_field.jl` 中添加接口**
添加以下接口并导出。建议同时补充文档。

```julia
function add_hex_anis(sim::AbstractSim, K1=0, K2=0, K3=0, name="hex")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    T = Float[]
    anis = HexagonalAnisotropy(T(K1), T(K2), T(K3), field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return anis
end
```

#### 5. **在 `src/micro/kernels.jl` 中实现 kernel**
```julia
@kernel function hexagonal_anisotropy_kernel!(@Const(m), h, energy, K1::T, K2::T, K3::T, 
                                            @Const(mu0_Ms), volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = mu0_Ms[id]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        Ms_inv::T = 1.0 / Ms_local
        @inbounds mx = m[j + 1]
        @inbounds my = m[j + 2]
        @inbounds mz = m[j + 3]
        @inbounds h[j + 1] = -6*K3*Ms_inv*(mx^5-10*mx^3*my^2+5*mx*my^4)
        @inbounds h[j + 2] = -6*K3*Ms_inv*(-5*mx^4*my+10*mx^2*my^3-my^5)
        @inbounds h[j + 3] = 2*mz*Ms_inv*(K1 + 2*K2*(1-mz*mz))
        @inbounds energy[id] = (K1*(1-mz*mz) + K2*(1-mz*mz)^2 + K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)) * volume
    end
end
```

#### 6. **在 `src/micro/field.jl` 中实现 `effective_field` 函数**
```julia
function effective_field(anis::HexagonalAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    hexagonal_anisotropy_kernel!(default_backend[])(spin, anis.field, anis.energy, anis.K1, 
                                            anis.K2, anis.K3, sim.mu0_Ms, volume; ndrange=N)

    return nothing
end
```

#### 7. **添加单元测试（`test/test_anis.jl`）**
在 `test_hex_anis` 函数中，对给定磁化计算有效场，并与自动微分得到的结果比较。

```julia
function hexagonal_energy(m, K1, K2, K3)
    mx, my, mz = m[1], m[2], m[3]
    return K1*(1-mz*mz) + K2*(1-mz*mz)^2 + K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)
end

function test_hex_anis()
    mesh = FDMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)
    m0 = (0.7, -0.4, 1.2)
    init_m0(sim, m0; norm=false)

    K1, K2, K3 = 1.23e2, 3.7e3, 6.9e2
    anis = add_hex_anis(sim, K1=K1, K2=K2, K3=K3)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    dEmx, dEmy, dEmz = hexagonal_energy_gradient(m0, K1, K2, K3)
    expected = [ -dEmx, -dEmy, -dEmz ] ./ (MicroMagnetic.mu_0*Ms)

    @test isapprox(field[1:3], expected)
    @test isapprox(energy[1]*1e27, hexagonal_energy(m0, K1, K2, K3), rtol=1e-5)
end
```

#### 8. **补充文档**
建议为新的实现补充文档，方便其他用户理解和使用。docstring 必须真正挂载,以下三条都踩过坑:

- 收尾的 `"""` 与 `function` 定义之间**不能有空行**——空行会静默断开挂载(在 Julia 1.12 上已验证);
- 用普通 `"""..."""` 或 **`@doc raw"""..."""`** 编写含公式的 docstring;裸 `raw"""..."""` 只是普通字符串,永远不会挂载;
- docstring 挂载到**紧随其后**的表达式——若在 docstring 和目标函数之间插入了辅助函数,docstring 会挂到辅助函数头上。

写完后用 `DOCS_DRAFT=true julia --project=docs docs/make.jl` 构建一次,检查有没有
`no docs found for ...` 警告,并把条目加进 `docs/src/api.md`。


### 关键命令速查
以下是本指南用到的关键 Git 命令速查：

- Fork：（通过 GitHub 网页）
- 克隆：`git clone https://github.com/yourusername/MicroMagnetic.jl`
- （可选）设置 upstream：`git remote add upstream https://github.com/magneticsimulation/MicroMagnetic.jl`
- 新建分支：`git checkout -b my-new-feature`
- 暂存改动：`git add .`
- 提交：`git commit -m "message"`
- 推送：`git push origin my-new-feature`
- （可选）拉取更新：`git fetch upstream`
- （可选）合并更新：`git merge upstream/main`
- 激活本地开发：`Pkg.develop(path="/path/to/MicroMagnetic.jl")`

欢迎探索、实验，一起把 **MicroMagnetic.jl** 做得更好！
