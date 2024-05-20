module oneAPIExt

import MicroMag
using oneAPI

function set_oneApi_backend()
    MicroMag.all_backends[3] = oneAPI.oneAPIBackend()
    MicroMag.set_backend("intel")
    return nothing
end

function __init__()
    set_oneApi_backend()
end

end
