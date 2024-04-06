module AMDGPUExt

import JuMag
using AMDGPU

function set_amd_backend()
    JuMag.all_backends[2] = AMDGPU.ROCBackend()
    JuMag.set_backend("amd")
    return nothing
end

function __init__()
    set_amd_backend()
end

end
