using PackageCompiler

# Runs inside a GPU-enabled container (see docker/build_cuda.sh); bakes CUDA
# kernels for the GPU architecture visible at build time into
# /opt/micromagnetic.so. Uses the default depot environment, which is where
# the Dockerfile installs the packages.
create_sysimage(["MicroMagnetic", "CUDA", "CUDSS", "CairoMakie", "NPZ"];
    sysimage_path = "/opt/micromagnetic.so",
    precompile_execution_file = "/opt/gpu_workload.jl")

@info "sysimage written" path="/opt/micromagnetic.so" size=Base.stat("/opt/micromagnetic.so").size
