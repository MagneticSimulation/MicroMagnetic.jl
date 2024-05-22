module oneAPIExt

import MicroMagnetic
using oneAPI

function set_oneApi_backend()
    MicroMagnetic.all_backends[3] = oneAPI.oneAPIBackend()
    MicroMagnetic.set_backend("intel")
    return nothing
end

function __init__()
    set_oneApi_backend()
end

end
