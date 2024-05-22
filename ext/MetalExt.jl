module MetalExt

import MicroMagnetic
using Metal

function set_metal_backend()
    MicroMagnetic.all_backends[4] = Metal.MetalBackend()
    MicroMagnetic.set_backend("apple")
    return nothing
end

function __init__()
    set_metal_backend()
end

end
