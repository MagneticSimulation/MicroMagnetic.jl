module CUDAExt

using MicroMagnetic
using CUDA
using CUDA.CUSPARSE

CUDA.allowscalar(false)

function set_cuda_backend()
    MicroMagnetic.all_backends[1] = CUDA.CUDABackend()
    MicroMagnetic.set_backend("cuda")
    MicroMagnetic.DefaultSparseMatrixCSC[] = CuSparseMatrixCSC
    return nothing
end

function __init__()
    return set_cuda_backend()
end

end
