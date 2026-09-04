# 周期性边界条件

MicroMagnetic 在三个坐标轴的任意组合上支持周期性边界条件(PBC)。周期性只在
**网格上声明一次**,所有相互作用项自动遵守:交换与 DMI 通过周期邻居表,退磁场
通过专用周期求解器。这消除了"交换周期而退磁开放"的静默不一致。

## 声明周期性

```julia
mesh = FDMesh(; nx=64, ny=64, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9, pbc="xy")
```

`pbc` 接受 `"x"`、`"y"`、`"z"` 字符的任意组合(默认 `"open"`,开放边界)。
其余搭建流程不变:

```julia
sim = Sim(mesh; driver="LLG", name="skx_lattice")
set_Ms(sim, 8e5)
add_exch(sim, 1.3e-11)
add_dmi(sim, 3e-3; type="interfacial")
add_demag(sim)                 # 依据 mesh.pbc 自动选择求解器
```

## 退磁求解器

退磁求解器依据网格周期性自动选择:

| 网格周期性 | 求解器 | 方法 | 文献 |
| :--- | :--- | :--- | :--- |
| 开放 | open demag | padded-FFT(或 `fft=false` 时直接求和) | Newell |
| 单轴(`x`/`y`/`z`) | pbc1d | 图像和 + 解析尾部 + DC 列修复 | Lebecki, *J. Phys. D* **41**, 175005 (2008) |
| 双轴(`xy`/`xz`/`yz`) | pbc2d | 图像和 + 解析尾部 + DC 列修复 | Wang et al., *Comp. Mater. Sci.* **49**, 84 (2010) |
| 三轴(`xyz`) | pbc3d | 谱投影器(tin-foil 约定) | Mansuripur (1989) |

周期化核为图像和 `N^PBC(d) = Σ_images N(d + p·L)`。单轴与双轴下该级数绝对收敛,
"显式图像 + 闭式连续尾部"是严格意义的分解。三轴下级数只是条件收敛,pbc3d 改用
标准谱投影器(均匀模不携带场)。

### 图像数 `Ic` / `Jc`

pbc1d/pbc2d 的每周期轴显式图像数默认 `Ic = Jc = 4`,加大可提升精度、增加一次性
初始化开销:

```julia
add_demag(sim; Ic=8, Jc=8)   # 更多图像,更小的尾部误差
```

相对暴力图像和的典型场误差:pbc1d `Ic=2` → 6×10⁻⁵、`Ic=4` → 2.5×10⁻⁶;
pbc2d `Ic=2` → 1×10⁻⁴。

### 遗留的显式 macro 模式

遗留的截断图像 macro 模式仍然可用,给出其关键字时优先生效:

```julia
add_demag(sim; macroPBC=true, Nx=2, Ny=2)   # Nx/Ny/Nz = 重复数
add_demag(sim; Nx=4, Ny=4)                  # 无 macroPBC 时的 legacy 语义
```

!!! note
    `fft=false` 会退化为直接求解器 + 截断图像;真正的 PBC 需要 FFT 管线(默认)。

## 一致性保证

- 周期轴只存在于网格上;`add_exch`、`add_dmi`、`add_demag` 全部从 `sim.mesh`
  读取,没有需要单独设置的逐项周期性关键字。
- 显式退磁关键字与网格周期轴矛盾时(例如在开放轴上给 `Nx`)会触发警告;
  在**开放**网格上传入则是静默的——这正是隔离测试 legacy macro 求解器的
  受支持用法。
- 周期网格上手工扫描邻居时,请使用带回绕的 `indexpbc(i, j, k, nx, ny, nz,
  xperiodic, yperiodic)` 而非 `index`。

!!! info "实现说明"
    周期退磁求解器的数学细节、符号约定与验证方法论见英文开发者文档
    `dev/demag.md`(Demag Internals)。
