## 安装

只要安装了 Julia（<http://julialang.org/downloads/>），安装 MicroMagnetic 非常简单，Windows、Linux 和 Mac 上都一样容易。

在 [Julia](http://julialang.org) 中，可以使用 Julia 包管理器轻松安装软件包。在 Julia REPL 中输入 `]` 进入 Pkg 模式，然后运行：

```julia
pkg> add MicroMagnetic
```

或者等价地：

```julia
julia> using Pkg;
julia> Pkg.add("MicroMagnetic")
```

要安装最新的开发版本：

```julia
pkg> add MicroMagnetic#master
```

要启用 GPU 支持，需要安装以下软件包之一：

!!! note "GPU 支持"
    | GPU 厂商               | Julia 软件包                                       |
    | :--------------------: | :------------------------------------------------: |
    | NVIDIA                 | [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)     |
    | AMD                    | [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl) |
    | Intel                  | [oneAPI.jl](https://github.com/JuliaGPU/oneAPI.jl) |
    | Apple                  | [Metal.jl](https://github.com/JuliaGPU/Metal.jl)   |

例如，对于 NVIDIA GPU 可以安装 `CUDA`：

```julia
pkg> add CUDA
```

安装后，输入 `using MicroMagnetic` 时会看到类似如下的信息：

```
julia> using MicroMagnetic
julia> using CUDA
Precompiling CUDAExt
  1 dependency successfully precompiled in 8 seconds. 383 already precompiled.
[ Info: Switch the backend to CUDA.CUDAKernels.CUDABackend(false, false)
```


## 从 Python 调用 MicroMagnetic.jl

借助 [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl)，可以无缝地从 Python 运行 MicroMagnetic.jl。

下面是一个**标准问题 #4** 的示例配置：

```python
from juliacall import Main as jl
jl.seval("using MicroMagnetic")

# 定义模拟参数
args = {
    "name": "std4",
    "task_s": ["relax", "dynamics"],                   # 要执行的任务列表
    "mesh": jl.FDMesh(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9),  # Julia 的 FDMesh 对象
    "Ms": 8e5,                                         # 饱和磁化强度 (A/m)
    "A": 1.3e-11,                                      # 交换刚度常数 (J/m)
    "demag": True,                                     # 启用退磁场
    "m0": (1, 0.25, 0.1),                              # 初始磁化矢量
    "alpha": 0.02,                                     # Gilbert 阻尼系数
    "steps": 100,                                      # 动力学模拟步数
    "dt": 0.01 * jl.ns,                                # 时间步长 (0.01 ns)
    "stopping_dmdt": 0.01,                             # 弛豫的停止判据
    "dynamic_m_every": 1,                              # 每步保存磁化
    "H_s": [(0, 0, 0), (-24.6 * jl.mT, 4.3 * jl.mT, 0)]  # 外磁场序列
}

# 运行模拟
sim = jl.sim_with(**args)
```
