module CUDAExt

using MicroMagnetic: MicroMagnetic
using CUDA

CUDA.allowscalar(false)

function set_cuda_backend()
    MicroMagnetic.all_backends[1] = CUDA.CUDABackend()
    MicroMagnetic.set_backend("cuda")
    return nothing
end

function __init__()
    return set_cuda_backend()
end

end
