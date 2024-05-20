module oneAPIExt

import NuMag
using oneAPI

function set_oneApi_backend()
    NuMag.all_backends[3] = oneAPI.oneAPIBackend()
    NuMag.set_backend("intel")
    return nothing
end

function __init__()
    set_oneApi_backend()
end

end
