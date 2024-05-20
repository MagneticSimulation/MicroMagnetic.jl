module CUDAExt

import NuMag
using CUDA

CUDA.allowscalar(false)

function set_cuda_backend()
    NuMag.all_backends[1] = CUDA.CUDABackend()
    NuMag.set_backend("cuda")
    return nothing
end

function __init__()
    set_cuda_backend()
end

end
