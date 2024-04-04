module AMDGPUExt

import JuMag
using AMDGPU

function set_amd_backend()
    JuMag.backend[] = AMDGPU.ROCBackend()
    @info("AMDGPU backend is used!")
    return nothing
end

function __init__()
    set_amd_backend()
end

end
