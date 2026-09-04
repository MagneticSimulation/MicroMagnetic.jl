# 容器

本指南介绍如何通过 Docker 和 Singularity/Apptainer 使用 MicroMagnetic.jl 容器镜像。

## 镜像概览

提供两个镜像：

| 镜像 | 内容 | 压缩拉取体积 |
|---|---|---|
| `ghcr.io/magneticsimulation/micromagnetic.jl:latest` | CUDA 版：Julia 1.12、CUDA 13.3 运行时与编译器（CUDA.jl artifact）、MicroMagnetic、CairoMakie、CUDSS、NPZ、**内置 sysimage** | 3.1 GB |
| `ghcr.io/magneticsimulation/micromagnetic.jl:cpu` | CPU 版：Julia 1.12、MicroMagnetic、CairoMakie、NPZ、**内置 sysimage** | 1.4 GB |

两个镜像都用 [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl) 构建了系统镜像
（`/opt/micromagnetic.so`），其中预编译了 MicroMagnetic、其依赖以及常用的模拟内核。
冷容器无需任何 depot 预热：在 64 核 A100 工作站上实测，CUDA 镜像 26 秒即可完成一条完整链路
（含 demag 的 GPU 模拟、绘图、导出视频），CPU 镜像 1.9 秒跑完完整 workload，全程零下载。

两个镜像的入口点均为 `julia -J /opt/micromagnetic.so`：脚本直接按路径传入（见下文示例），
不带参数的 `docker run -it` 会启动一个 sysimage 生效的 Julia REPL。

### 可用标签

| 标签 | 内容 |
|---|---|
| `:latest` | CUDA 版（默认） |
| `:cuda` | `:latest` 的别名 |
| `:cpu` | 纯 CPU 版 |
| `:0.7.0`、`:0.7.0-cpu` | 冻结在 MicroMagnetic v0.7.0 的镜像 |

CUDA 镜像重建时会自动打上 `:cuda` 以及 `Project.toml` 中包版本对应的标签
（见 `docker/build_cuda.sh`）。

说明：

- CUDA 运行时与编译器库来自 CUDA.jl 的 artifact 分发；NVIDIA 驱动在运行时由
  nvidia-container-toolkit 注入，**不会**打进镜像。
- 宿主机驱动需支持 CUDA 13（驱动 ≥ 580）。
- sysimage 中的 CUDA kernel 是针对**构建镜像时所用 GPU 的计算能力**预编译的
  （例如 A100 为 sm_80）。在其他代际的 GPU 上，受影响的 kernel 会在首次使用时
  JIT 编译一次。

## 构建镜像

CUDA sysimage 必须在装有 NVIDIA GPU 和 nvidia-container-toolkit 的机器上构建，
因为 `docker build` 无法访问 GPU。
[`docker/build_cuda.sh`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/docker/build_cuda.sh)
把三个步骤串了起来：`docker build`（依赖）→ `docker run --gpus all`
（烘 sysimage + 全路径探针）→ `docker commit`（成品镜像，自动打 `:latest`、`:cuda`
以及 `Project.toml` 版本号标签）：

```bash
./docker/build_cuda.sh ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

CPU 镜像无需 GPU，一步构建：

```bash
docker build -f docker/Dockerfile.cpu -t ghcr.io/magneticsimulation/micromagnetic.jl:cpu .
```

## Docker 使用

### 前置条件

- [Docker](https://docs.docker.com/engine/install/)
- GPU 支持需要 [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
  以及驱动 ≥ 580

### 拉取镜像

```bash
docker pull ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

国内用户可通过 NJU 镜像源拉取：

```bash
docker pull ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
```

### 运行模拟脚本

入口点为 `julia -J /opt/micromagnetic.so`，直接传脚本路径即可（**不要**再写前导 `julia`）：
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  /workspace/run_simulation.jl
```

纯 CPU 主机使用 `:cpu` 镜像并去掉 `--gpus all`：
```bash
docker run --rm -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:cpu \
  /workspace/run_simulation.jl
```

### 交互式会话

```bash
docker run -it --rm --gpus all ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

附加 Julia 参数写在镜像名之后、脚本之前：
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  -t 8 /workspace/run_simulation.jl
```

### 不带 sysimage 的原生 Julia

调试包问题时，可用原生系统镜像启动 Julia：
```bash
docker run -it --rm --gpus all --entrypoint julia \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

### 安装额外的包

镜像内置的 depot（`/usr/local/share/julia`）已包含运行所需的一切。若要永久添加包，
在容器里安装后提交：

```bash
docker run -it --name mm-extra ghcr.io/magneticsimulation/micromagnetic.jl:latest
# 容器内：julia -e 'using Pkg; Pkg.add("MyPackage")'，然后 exit
docker commit mm-extra my-micromagnetic && docker rm mm-extra
```

或派生一个镜像：

```dockerfile
FROM ghcr.io/magneticsimulation/micromagnetic.jl:latest
RUN julia -e 'using Pkg; Pkg.add("MyPackage")'
```

## Singularity / Apptainer 使用

Singularity/Apptainer 是 HPC 集群上的标准容器运行时，无需 root 权限，并与 SLURM 集成。

`--writable-tmpfs` 是 GPU 运行的必需参数：CUDA 栈启动时会向（只读的）镜像 depot 写入
scratch 使用日志，该参数提供一层临时的可写层。下述 GPU 链路（含 demag 的模拟、绘图、
视频导出）已用这条命令实测通过。

### 拉取镜像

```bash
singularity pull micromagnetic.sif docker://ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

### 运行模拟脚本

Apptainer 的 `exec` 会用给定命令**替代**镜像入口点，因此需要用 `-J` 显式指定 sysimage：
```bash
singularity exec \
  --writable-tmpfs \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  julia -J /opt/micromagnetic.so /workspace/run_simulation.jl
```

### 交互式会话

`run` 会遵循镜像入口点，启动 sysimage 生效的 REPL：
```bash
singularity run --nv --writable-tmpfs micromagnetic.sif
```

### SLURM 作业示例

```bash
#!/bin/bash
#SBATCH --gpus=1
#SBATCH --time=01:00:00

singularity exec \
  --writable-tmpfs \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  julia -J /opt/micromagnetic.so /workspace/run_simulation.jl
```

## 常见问题

### GPU 未被识别

- **Docker**：确认带了 `--gpus all`
- **Singularity**：确认带了 `--nv`
- 先在宿主机确认驱动正常（`nvidia-smi` 可用）

### CUDA 驱动过旧

CUDA 镜像内置 CUDA 13.3 运行时库，宿主机驱动需 ≥ 580。若要在更旧的驱动上运行，
可在 `docker/Dockerfile.cuda` 中重建镜像时用官方 API 锁定旧版运行时
（在包安装之后、sysimage 构建之前的任意位置）：

```dockerfile
RUN julia -e 'using CUDA; CUDA.set_runtime_version!(v"12.6")'
```

CUDA 12.6 的 artifact 可运行于驱动 ≥ 525。

### 在不同 GPU 上 kernel 重新编译

sysimage 中的 CUDA kernel 是按构建时 GPU 的计算能力烘焙的。在其他代际的 GPU 上，
kernel 会 JIT 编译一次（并缓存进 depot）；如需消除，可在装该 GPU 的机器上重建镜像。

### HPC 集群上找不到 singularity

- 先试 `module load singularity`（或 `module load apptainer`）
- 否则按 https://apptainer.org/ 的说明安装
