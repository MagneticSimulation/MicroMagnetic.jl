# Container

This guide covers using the MicroMagnetic.jl containers with Docker and Singularity/Apptainer.

## Image Overview

Two images are provided:

| Image | Contents | Compressed pull size |
|---|---|---|
| `ghcr.io/magneticsimulation/micromagnetic.jl:latest` | CUDA build: Julia 1.12, CUDA 13.3 runtime & compiler (CUDA.jl artifacts), MicroMagnetic, CairoMakie, CUDSS, NPZ, baked sysimage | 3.1 GB |
| `ghcr.io/magneticsimulation/micromagnetic.jl:cpu` | CPU build: Julia 1.12, MicroMagnetic, CairoMakie, NPZ, baked sysimage | 1.4 GB |

Both images contain a [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl)
system image at `/opt/micromagnetic.so` with MicroMagnetic, its dependencies and the hot
simulation kernels pre-compiled. No Julia depot warm-up is needed: on a 64-core A100
workstation a cold container runs a full path (GPU simulation with demag, plotting, movie
export) in 26 s, and the CPU image runs its full workload in 1.9 s — with zero downloads.

The entrypoint of both images is `julia -J /opt/micromagnetic.so`: scripts are passed by
path (see the examples below), and a bare `docker run -it` starts a Julia REPL with the
sysimage active.

### Available tags

| Tag | Contents |
|---|---|
| `:latest` | CUDA build (default) |
| `:cuda` | alias of `:latest` |
| `:cpu` | CPU-only build |
| `:0.7.0`, `:0.7.0-cpu` | images frozen at MicroMagnetic v0.7.0 |

Rebuilds of the CUDA image are automatically tagged `:cuda` and with the package version from
`Project.toml` (see `docker/build_cuda.sh`).

Notes:

- CUDA runtime and compiler libraries come from the CUDA.jl artifact distribution; the NVIDIA
  driver is injected at runtime by nvidia-container-toolkit and is never baked into the image.
- The host must have an NVIDIA driver that supports CUDA 13 (driver ≥ 580).
- CUDA kernels are baked into the sysimage for the **compute capability of the GPU present at
  build time** (e.g. sm_80 for A100). On other GPU generations the affected kernels are
  JIT-compiled once on first use.

## Docker Usage

### Prerequisites

- [Docker](https://docs.docker.com/engine/install/)
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
  for GPU support, and an NVIDIA driver ≥ 580

### Pull the image

```bash
docker pull ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

Users in China can pull through the NJU mirror:
```bash
docker pull ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
```

### Run a simulation script

The entrypoint is `julia -J /opt/micromagnetic.so`, so pass the script path directly (no
leading `julia`):
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  /workspace/run_simulation.jl
```

On CPU-only hosts use the `:cpu` image and drop `--gpus all`:
```bash
docker run --rm -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:cpu \
  /workspace/run_simulation.jl
```

### Interactive session

```bash
docker run -it --rm --gpus all ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

Extra Julia options go after the image name, before the script:
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  -t 8 /workspace/run_simulation.jl
```

### Plain Julia without the sysimage

To start Julia with the stock system image (e.g. for debugging package issues):
```bash
docker run -it --rm --gpus all --entrypoint julia \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

### Installing extra packages

The depot baked into the image (`/usr/local/share/julia`) already contains everything needed
to run. To add a package permanently, commit a container after installing it:
```bash
docker run -it --name mm-extra ghcr.io/magneticsimulation/micromagnetic.jl:latest
# inside the container:
#   julia -e 'using Pkg; Pkg.add("MyPackage")'
#   exit
docker commit mm-extra my-micromagnetic && docker rm mm-extra
```
or derive an image:
```dockerfile
FROM ghcr.io/magneticsimulation/micromagnetic.jl:latest
RUN julia -e 'using Pkg; Pkg.add("MyPackage")'
```

## Singularity / Apptainer Usage

Singularity/Apptainer is the standard container runtime on HPC clusters. It runs without
root privileges and integrates with SLURM.

`--writable-tmpfs` is required for GPU runs: the CUDA stack writes a scratch usage log into
the (read-only) image depot at startup, and this flag provides an ephemeral writable layer.
The GPU path (simulation with demag, plotting, movie export) is verified with this exact
invocation.

### Pull the image

```bash
singularity pull micromagnetic.sif docker://ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

### Run a simulation script

Apptainer's `exec` runs the given command *instead of* the image entrypoint, so the sysimage
must be passed explicitly with `-J`:
```bash
singularity exec \
  --writable-tmpfs \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  julia -J /opt/micromagnetic.so /workspace/run_simulation.jl
```

### Interactive session

`run` follows the image entrypoint and starts the REPL with the sysimage:
```bash
singularity run --nv --writable-tmpfs micromagnetic.sif
```

### SLURM Job Example

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

## Building the Images

The CUDA sysimage must be built on a machine with an NVIDIA GPU and
`nvidia-container-toolkit`, because `docker build` cannot access GPUs.
[`docker/build_cuda.sh`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/docker/build_cuda.sh)
does the three steps for you: `docker build` (dependencies) → `docker run --gpus all`
(sysimage + full-path probe) → `docker commit` (final image, tagged `:latest`, `:cuda` and
with the version from `Project.toml`):

```bash
./docker/build_cuda.sh ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

The CPU image needs no GPU and builds in one step:

```bash
docker build -f docker/Dockerfile.cpu -t ghcr.io/magneticsimulation/micromagnetic.jl:cpu .
```

## Common Issues

### GPU not detected

- **Docker**: make sure `--gpus all` is present
- **Singularity**: make sure `--nv` is present
- Verify that NVIDIA drivers are installed on the host system (`nvidia-smi` works)

### CUDA driver too old

The CUDA image ships CUDA 13.3 runtime libraries, so the host driver must be ≥ 580. To run on
older drivers, rebuild the image with an older CUDA runtime: pin the version in
`docker/Dockerfile.cuda` with the official API (any time after the packages are installed and
before the sysimage build),

```dockerfile
RUN julia -e 'using CUDA; CUDA.set_runtime_version!(v"12.6")'
```

CUDA 12.6 artifacts run on drivers ≥ 525.

### Kernels recompile on a different GPU

The sysimage bakes CUDA kernels for the build-time GPU's compute capability. On a different
GPU generation the kernels JIT-compile once (and are cached in the depot); rebuild the image
on a machine with that GPU to bake them.

### singularity not found on HPC clusters

- Try `module load singularity` (or `module load apptainer`)
- Otherwise install it from https://apptainer.org/
