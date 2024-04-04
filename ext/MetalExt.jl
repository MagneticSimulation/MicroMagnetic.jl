module MetalExt

import JuMag
using Metal

function set_metal_backend()
    JuMag.backend[] = Metal.MetalBackend()
    @info("Apple Metal backend is used!")
    return nothing
end

function __init__()
    set_metal_backend()
end

end
