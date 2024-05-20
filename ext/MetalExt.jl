module MetalExt

import NuMag
using Metal

function set_metal_backend()
    NuMag.all_backends[4] = Metal.MetalBackend()
    NuMag.set_backend("apple")
    return nothing
end

function __init__()
    set_metal_backend()
end

end
