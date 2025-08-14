module CUDAExt

using MicroMagnetic
using CUDA
using CUDA.CUSPARSE

CUDA.allowscalar(false)

function set_cuda_backend()
    MicroMagnetic.all_backends[1] = CUDA.CUDABackend()
    MicroMagnetic.set_backend("cuda")
    MicroMagnetic.GPUSparseMatrixCSC[] = CuSparseMatrixCSC
    MicroMagnetic.GPUSparseMatrixCSR[] = CuSparseMatrixCSR
    MicroMagnetic.GPUMatrix[] = CuMatrix
    return nothing
end

function __init__()
    return set_cuda_backend()
end

end
