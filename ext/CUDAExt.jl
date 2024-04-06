module CUDAExt

import JuMag
using CUDA

CUDA.allowscalar(false)

function set_cuda_backend()
    JuMag.all_backends[1] = CUDA.CUDABackend()
    JuMag.set_backend("cuda")
    return nothing
end

function __init__()
    set_cuda_backend()
end

end
