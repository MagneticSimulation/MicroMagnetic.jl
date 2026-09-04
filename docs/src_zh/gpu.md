# GPU 与精度

MicroMagnetic.jl 默认在 CPU 上运行,并通过
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) 支持
NVIDIA、AMD、Intel 和 Apple GPU。本页说明三个会话级设置——**后端**、**精度**与
**线程组大小**——以及一条需要牢记的规则:

!!! warning "设置是构造期默认值"
    `set_backend`、`set_precision` 和 `set_groupsize` 只影响**之后**创建的模拟。
    `Sim` 在构造时把精度冻结进类型、把后端冻结进数据;之后再改全局设置,
    不会触碰已经存在的模拟。

## 后端

按[安装页](install.md)的说明装好 GPU 包后,用 `@using_gpu` 宏加载——它会自动
导入当前环境中可用的后端包:

```julia
using MicroMagnetic
@using_gpu()   # 有 CUDA/AMDGPU/oneAPI/Metal 则加载之,否则留在 CPU
```

该宏会自动切换默认后端。如需手动控制:

```julia
set_backend("cuda")   # 之后的分配进入 GPU
set_backend("cpu")    # 回到 CPU
```

已有模拟存在时调用 `set_backend` 会打印警告——正在运行的模拟数据保持原地不动。

## 精度

默认工作精度是 `Float64`。单精度在消费级 GPU 上通常能翻倍吞吐
(GeForce 卡 FP64 吞吐只有 FP32 的 1/32),并把显存占用减半:

```julia
set_precision(Float32)   # 这行之后创建的模拟使用 Float32
sim = Sim(mesh; driver="LLG")
```

允许的取值为 `Float64`、`Float32` 和 `AbstractFloat`:

- `Float64` / `Float32` —— 生产运行的具象精度;
- `AbstractFloat` —— **符号模式**,只有本征值分析工作流使用,允许参数携带符号项。
  常规模拟在该模式下会失败,尤其在 GPU 上(设备数组不接受抽象元素类型)。

!!! note "先设精度,再建模拟"
    精度在构造时冻结进模拟类型,所以要在 `Sim(...)` / `CubicMesh(...)` /
    `FDMesh(...)` **之前**调用 `set_precision`。不支持同一会话内混用精度。

### 本征值分析的例外

本征值分析需要符号模式,它会修改**整个会话**的全局精度且不会自动恢复。
做完本征值计算后,在创建后续模拟之前先恢复具象精度:

```julia
# ... set_precision(AbstractFloat) 下的本征值分析 ...
set_precision(Float64)   # 下一个模拟之前恢复
```

完整工作流见[本征模分析](eigen/eigenmodes.md)页。

## 线程组大小

`set_groupsize` 控制 GPU kernel 发射时的线程组大小(默认 256)。它是纯性能
旋钮——不影响任何计算结果:

```julia
set_groupsize(128)   # 只有 profiling 建议时才调
```

## 查看当前状态

`using MicroMagnetic` 后的启动信息会报告当前后端。脚本里可以通过模拟本身
确认冻结的精度:

```julia
sim = Sim(mesh; driver="LLG")
eltype(sim.spin)                    # 该模拟的工作精度
```

!!! tip "先测量,再调参"
    任何性能相关的事——组大小、精度、后端——都先用 `MicroMagnetic.timer`
    前后测量;报告长什么样见[基础页](basics.md)的 Timings 一节。
