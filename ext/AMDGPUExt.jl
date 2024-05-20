module AMDGPUExt

import NuMag
using AMDGPU

function set_amd_backend()
    NuMag.all_backends[2] = AMDGPU.ROCBackend()
    NuMag.set_backend("amd")
    return nothing
end

function __init__()
    set_amd_backend()
end

end
