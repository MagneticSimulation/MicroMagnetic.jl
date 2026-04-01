# Container

This guide covers using the MicroMagnetic.jl container with Docker and Singularity/Apptainer.

## Image Overview

The [container image](https://github.com/MagneticSimulation/MicroMagnetic.jl/pkgs/container/micromagnetic.jl) includes:

- Julia 1.10 with CUDA.jl pre-installed
- MicroMagnetic.jl (development version from `master`)
- CairoMakie (simplied version) for visualization

**Image size**: ~4.2 GB

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
<<<<<<< HEAD
docker pull ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
=======
docker pull  docker://ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
>>>>>>> 33832f35f6700d5947e12a8dd596d2a3003a38cb
```

**Run a simulation script:**
```bash
docker run --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  julia /workspace/run_simulation.jl
```

**Start an interactive Julia session:**
```bash
docker run -it --rm --gpus all -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

**Persistent cache (recommended for faster repeated runs):**

To reuse precompiled packages across multiple runs, bind a persistent directory to `/depot`:

```bash
mkdir -p ~/julia_depot

docker run --rm --gpus all \
  -v ~/julia_depot:/depot \
  -v $(pwd):/workspace \
  ghcr.io/magneticsimulation/micromagnetic.jl:latest \
  julia /workspace/run_simulation.jl
```

Without this, precompiled packages are discarded after each run, causing slower startup.

## Singularity Usage

Singularity/Apptainer is the standard container runtime on HPC clusters. It runs without root privileges and integrates with SLURM.

**Important**: Singularity mounts containers as read-only by default. You **must** bind a writable directory to `/depot` for Julia to function.

### Prerequisites

- SingularityCE 3.0+ or Apptainer 1.0+
- GPU support requires the `--nv` flag

### Commands

**Pull the image:**
```bash
singularity pull micromagnetic.sif docker://ghcr.io/magneticsimulation/micromagnetic.jl:latest
```

Alternatively, you can pull the image from `ghcr.nju.edu.cn` to speed up the download:
```bash
singularity pull micromagnetic.sif docker://ghcr.nju.edu.cn/magneticsimulation/micromagnetic.jl:latest
```

**Run a script (with required depot binding):**
```bash
mkdir -p /work/${USER}/julia_depot

singularity exec \
  --bind /work/${USER}/julia_depot:/depot \
  --nv \
  micromagnetic.sif \
  julia /path/to/script.jl
```

**Interactive session:**
```bash
singularity shell \
  --bind /work/${USER}/julia_depot:/depot \
  --nv \
  micromagnetic.sif
```

### SLURM Job Example

```bash
#!/bin/bash
#SBATCH --gpus=1
#SBATCH --time=01:00:00

export JULIA_DEPOT=/work/${USER}/julia_depot
mkdir -p $JULIA_DEPOT

singularity exec \
  --bind $JULIA_DEPOT:/depot \
  --bind $(pwd):/workspace \
  --nv \
  micromagnetic.sif \
  julia /workspace/run_simulation.jl
```
---

## Common Issues

### Singularity: read-only file system

**Error**: `read-only file system` when running Julia.

**Solution**: Bind a writable directory to `/depot`:
```bash
mkdir -p ~/julia_depot
singularity exec --bind ~/julia_depot:/depot ...
```

### GPU not detected

- **Docker**: Ensure `--gpus all` is included
- **Singularity**: Ensure `--nv` is included
- Verify NVIDIA drivers are installed on the host system

### singularity not found in HPC clusters
- try `module load singularity` to load the Singularity module
- install it by following the instructions on the https://apptainer.org/
