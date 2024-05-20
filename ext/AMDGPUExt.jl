module AMDGPUExt

import MicroMag
using AMDGPU

function set_amd_backend()
    MicroMag.all_backends[2] = AMDGPU.ROCBackend()
    MicroMag.set_backend("amd")
    return nothing
end

function __init__()
    set_amd_backend()
end

end
