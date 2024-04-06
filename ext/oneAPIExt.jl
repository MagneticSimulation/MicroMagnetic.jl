module oneAPIExt

import JuMag
using oneAPI

function set_oneApi_backend()
    JuMag.all_backends[3] = oneAPI.oneAPIBackend()
    JuMag.set_backend("intel")
    return nothing
end

function __init__()
    set_oneApi_backend()
end

end
