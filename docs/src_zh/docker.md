# 容器

本指南介绍如何通过 Docker 和 Singularity/Apptainer 使用 MicroMagnetic.jl 容器镜像。

## 镜像概览

提供两个镜像：

| 镜像 | 内容 | 典型大小 |
|---|---|---|
| `micromagnetic.jl:latest`（CUDA） | Julia + CUDA toolkit（系统） + MicroMagnetic + CairoMakie，**内置 sysimage** | ~4 GB |
| `micromagnetic.jl:cpu` | Julia + MicroMagnetic + CairoMakie，**内置 sysimage** | ~1.5 GB |

两个镜像都带有用 [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl) 构建的
系统镜像（`/opt/micromagnetic.so`），其中预编译了 MicroMagnetic、其依赖以及常用的模拟内核。
启动不再依赖预热好的 Julia depot：`using MicroMagnetic` 只需几秒而不是约 30–60 秒，
首次 GPU kernel 启动也是即时的，无需为每个 kernel 支付约 9 秒的 JIT 编译
（在 A100、Julia 1.12 上测得；见下文）。

注意：

- sysimage 中的 CUDA kernel 是针对**构建镜像时所用 GPU 的计算能力**预编译的
  （例如 A100 为 sm_80）。在其他代际的 GPU 上，受影响的 kernel 会在首次使用时
  像以前一样 JIT 编译一次。
- CUDA 镜像设置了 `local_toolkit` preference，因此 CUDA.jl 使用基础镜像中安装的
  toolkit，而不是把数 GB 的运行时 artifact 下载进 Julia depot。找不到系统 toolkit 时
  会自动回退到 artifact。
- 宿主机需要安装与镜像 CUDA toolkit 兼容的 NVIDIA 驱动
  （CUDA 12.6 ⇔ 驱动 ≥ 560）。

## 构建镜像

CUDA sysimage 必须在有 NVIDIA GPU 的机器上构建（并安装 `nvidia-container-toolkit`），
因为 `docker build` 无法访问 GPU。[docker/build_cuda.sh](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/docker/build_cuda.sh)
帮你完成三步：`docker build`（依赖）→ `docker run --gpus all`
（sysimage）→ `docker commit`（最终镜像）：

```bash
./docker/build_cuda.sh ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

CPU 镜像一步即可构建（无需 GPU）：

```bash
docker build -f docker/Dockerfile.cpu -t ghcr.io/magneticsimulation/micromagnetic.jl:cpu .
```

## Docker 用法

### 前置条件

- [Docker Desktop](https://www.docker.com/products/docker-desktop/)
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)（GPU 支持需要）

### 常用命令

**拉取镜像：**
```bash
docker pull ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

也可以从 `ghcr.nju.edu.cn` 拉取镜像以加速下载：
```bash
docker pull ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
```

**运行模拟脚本。** 入口点是 `julia -J /opt/micromagnetic.so`，直接传脚本路径即可
（不要带前导 `julia`）：
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  /workspace/run_simulation.jl
```

**启动交互式 Julia 会话（使用 sysimage）：**
```bash
docker run -it --rm --gpus all ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**附加 Julia 选项**（线程、project 等）可以同样追加：
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  -t 8 /workspace/run_simulation.jl
```

**不使用 sysimage 的纯 Julia**（例如要 `Pkg.add` 到 depot 里）：
```bash
docker run -it --rm --gpus all --entrypoint julia \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**持久化 depot（可选）。** 有了内置 sysimage，即使 depot 是全新的，启动也很快，
因此不再需要持久化 depot。如果想安装能在多次运行之间保留的额外软件包，
可以挂载一个目录并把 `JULIA_DEPOT_PATH` 指向它：
```bash
mkdir -p ~/julia_depot
docker run --rm --gpus all \
  -e JULIA_DEPOT_PATH=/depot \
  -v ~/julia_depot:/depot \
  -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  /workspace/run_simulation.jl
```

## Singularity 用法

Singularity/Apptainer 是 HPC 集群上的标准容器运行时。它无需 root 权限即可运行，
并与 SLURM 集成。

**重要**：Singularity 默认以只读方式挂载容器。如果 Julia 需要写入
（临时输出、Pkg 操作等），请绑定一个可写目录；内置 sysimage 本身不需要可写 depot。

### 前置条件

- SingularityCE 3.0+ 或 Apptainer 1.0+
- GPU 支持需要 `--nv` 选项

### 常用命令

**拉取镜像：**
```bash
singularity pull micromagnetic.sif docker://ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**运行脚本**（入口点已带 sysimage，直接传脚本路径）：
```bash
singularity exec \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  /workspace/run_simulation.jl
```

**交互式会话：**
```bash
singularity shell --nv micromagnetic.sif
```

### SLURM 作业示例

```bash
#!/bin/bash
#SBATCH --gpus=1
#SBATCH --time=01:00:00

singularity exec \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  /workspace/run_simulation.jl
```

---

## 常见问题

### 检测不到 GPU

- **Docker**：确保命令中包含 `--gpus all`
- **Singularity**：确保命令中包含 `--nv`
- 确认宿主机已安装 NVIDIA 驱动

### CUDA 驱动过旧

CUDA 镜像自带 CUDA 12.6 toolkit，宿主机驱动必须 ≥ 560。若驱动较旧，
可通过 `--build-arg CUDA_BASE=nvidia/cuda:12.4.1-devel-ubuntu22.04` 用旧版 toolkit 构建。

### 换 GPU 后 kernel 重新编译

sysimage 为构建时 GPU 的计算能力预编译了 CUDA kernel。在其他代际的 GPU 上，
kernel 会 JIT 编译一次（并缓存到 depot 中）；要在该 GPU 上预编译，请在装有该 GPU 的
机器上重新构建镜像。

### HPC 集群上找不到 singularity
- 尝试 `module load singularity` 加载 Singularity 模块
- 或按照 https://apptainer.org/ 的说明自行安装
