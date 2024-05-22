module AMDGPUExt

import MicroMagnetic
using AMDGPU

function set_amd_backend()
    MicroMagnetic.all_backends[2] = AMDGPU.ROCBackend()
    MicroMagnetic.set_backend("amd")
    return nothing
end

function __init__()
    set_amd_backend()
end

end
