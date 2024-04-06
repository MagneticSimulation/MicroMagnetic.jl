module MetalExt

import JuMag
using Metal

function set_metal_backend()
    JuMag.all_backends[4] = Metal.MetalBackend()
    JuMag.set_backend("apple")
    return nothing
end

function __init__()
    set_metal_backend()
end

end
