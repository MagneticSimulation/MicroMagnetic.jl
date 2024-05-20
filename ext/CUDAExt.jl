module CUDAExt

import MicroMag
using CUDA

CUDA.allowscalar(false)

function set_cuda_backend()
    MicroMag.all_backends[1] = CUDA.CUDABackend()
    MicroMag.set_backend("cuda")
    return nothing
end

function __init__()
    set_cuda_backend()
end

end
