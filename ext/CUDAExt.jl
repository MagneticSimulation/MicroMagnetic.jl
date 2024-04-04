module CUDAExt

import JuMag
using CUDA

CUDA.allowscalar(false)

function set_cuda_backend()
    JuMag.backend[] = CUDA.CUDABackend()
    @info("NVIDIA CUDA backend is used!")
    return nothing
end

function __init__()
    set_cuda_backend()
end

end
