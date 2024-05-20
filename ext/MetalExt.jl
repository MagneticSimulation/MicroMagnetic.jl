module MetalExt

import MicroMag
using Metal

function set_metal_backend()
    MicroMag.all_backends[4] = Metal.MetalBackend()
    MicroMag.set_backend("apple")
    return nothing
end

function __init__()
    set_metal_backend()
end

end
