# Container

This guide covers using the MicroMagnetic.jl containers with Docker and Singularity/Apptainer.

## Image Overview

Two images are provided:

| Image | Contents | Size (uncompressed) |
|---|---|---|
| `micromagnetic.jl:latest` (CUDA) | Julia 1.12 + CUDA.jl (runtime & compiler from artifacts) + MicroMagnetic + CairoMakie, **sysimage baked** | ~3.3 GB |
| `micromagnetic.jl:cpu` | Julia + MicroMagnetic + CairoMakie, **sysimage baked** | ~1.5 GB |

Cold-container full-path timings measured on a 64-core A100 workstation: CUDA
image 26 s (GPU simulation + plotting + movie export), CPU image 1.9 s
(atomistic + FDMesh simulations).

CUDA libraries come from CUDA.jl's artifact distribution and the NVIDIA driver is
injected by nvidia-container-toolkit at runtime, so the image carries no system
CUDA toolkit; the bulk is the depot artifacts plus the ~1.5 GB sysimage that
makes startup instant. For a small download, use the CPU image.

Both images ship a [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl) system
image (`/opt/micromagnetic.so`) with MicroMagnetic, its dependencies and the hot simulation
kernels pre-compiled. Startup no longer depends on a warmed Julia depot: `using MicroMagnetic`
takes a few seconds instead of ~30–60 s, and the first GPU kernel launch is fast instead of
paying a ~9 s JIT compile per kernel (measured on A100, Julia 1.12; see below).

Notes:

- CUDA kernels in the sysimage are baked for the **compute capability of the GPU present at
  image build time** (e.g. sm_80 for A100). On other GPU generations the affected kernels are
  JIT-compiled once on first use, as before.
- CUDA runtime & compiler libraries come from the CUDA.jl artifact distribution (CUDA 13.3);
  the NVIDIA driver is injected at runtime by nvidia-container-toolkit and never baked into
  the image.
- The host must have an NVIDIA driver that supports CUDA 13 (driver ≥ 580).

## Building the images

The CUDA sysimage must be built on a machine with an NVIDIA GPU (and `nvidia-container-toolkit`
installed), because `docker build` cannot access GPUs. [`docker/build_cuda.sh`](https://github.com/MagneticSimulation/MicroMagnetic.jl/blob/master/docker/build_cuda.sh)
does the three steps for you: `docker build` (dependencies) → `docker run --gpus all`
(sysimage + full-path probe) → `docker commit` (final image, tagged `:latest`, `:cuda` and
with the version from `Project.toml`):

```bash
./docker/build_cuda.sh ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

The CPU image builds in one step (no GPU needed):

```bash
docker build -f docker/Dockerfile.cpu -t ghcr.io/magneticsimulation/micromagnetic.jl:cpu .
```

### Available tags

| Tag | Contents |
|---|---|
| `:latest` | CUDA build (default) |
| `:cuda` | alias of `:latest` |
| `:cpu` | CPU-only build |
| `:0.7.0`, `:0.7.0-cpu` | images frozen at MicroMagnetic v0.7.0 |

Rebuilds of the CUDA image are automatically tagged `:cuda` and with the package version from
`Project.toml` (see `docker/build_cuda.sh`).

## Docker Usage

### Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop/)
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) for GPU support

### Commands

**Pull the image:**
```bash
docker pull ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

Alternatively, you can pull the image from `ghcr.nju.edu.cn` to speed up the download:
```bash
docker pull ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
```

**Run a simulation script.** The entrypoint is `julia -J /opt/micromagnetic.so`, so pass the
script path directly (no leading `julia`):
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  /workspace/run_simulation.jl
```

**CPU-only hosts:** use the `:cpu` image and drop `--gpus all`:
```bash
docker run --rm -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:cpu \
  /workspace/run_simulation.jl
```

**Start an interactive Julia session (with the sysimage):**
```bash
docker run -it --rm --gpus all ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**Extra Julia options** (threads, project, …) can be appended the same way:
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  -t 8 /workspace/run_simulation.jl
```

**Plain Julia without the sysimage** (e.g. to `Pkg.add` into the depot):
```bash
docker run -it --rm --gpus all --entrypoint julia \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**Installing extra packages.** The image's Julia depot (`/usr/local/share/julia`) already
contains everything needed to run — do **not** point `JULIA_DEPOT_PATH` at an empty mounted
directory: that hides the baked CUDA artifacts and registries and forces a multi-GB
re-download. To add packages permanently, commit a container after installing them,
```bash
docker run -it --name mm-extra ghcr.io/magneticsimulation/micromagnetic.jl:latest
# inside: julia -e 'using Pkg; Pkg.add("MyPackage")'; exit
docker commit mm-extra my-micromagnetic && docker rm mm-extra
```
or derive an image:
```dockerfile
FROM ghcr.io/magneticsimulation/micromagnetic.jl:latest
RUN julia -e 'using Pkg; Pkg.add("MyPackage")'
```

## Singularity Usage

Singularity/Apptainer is the standard container runtime on HPC clusters. It runs without root
privileges and integrates with SLURM.

**Important**: Singularity mounts containers as read-only by default. If Julia needs to write
(temporary output, Pkg operations), bind a writable directory; the baked sysimage itself needs
no writable depot.

### Prerequisites

- SingularityCE 3.0+ or Apptainer 1.0+
- GPU support requires the `--nv` flag

### Commands

**Pull the image:**
```bash
singularity pull micromagnetic.sif docker://ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**Run a script** (the entrypoint carries the sysimage, pass the script path directly):
```bash
singularity exec \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  /workspace/run_simulation.jl
```

**Interactive session:**
```bash
singularity shell --nv micromagnetic.sif
```

### SLURM Job Example

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

## Common Issues

### GPU not detected

- **Docker**: Ensure `--gpus all` is included
- **Singularity**: Ensure `--nv` is included
- Verify NVIDIA drivers are installed on the host system

### CUDA driver too old

The CUDA image uses CUDA 13.3 runtime libraries (from the CUDA.jl artifact); the
host driver must support CUDA 13 (≥ 580). For older drivers, build with an older
runtime via `--build-arg CUDA_BASE=nvidia/cuda:12.4.1-base-ubuntu22.04` and pin
the runtime version, or install an older CUDA.jl.

### Kernels recompile on a different GPU

The sysimage bakes CUDA kernels for the build-time GPU's compute capability. On a different GPU
generation the kernels JIT-compile once (and are cached in the depot); rebuild the image on a
machine with that GPU to bake them.

### singularity not found in HPC clusters
- try `module load singularity` to load the Singularity module
- install it by following the instructions on the https://apptainer.org/
