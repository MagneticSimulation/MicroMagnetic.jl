module oneAPIExt

import JuMag
using oneAPI

function set_oneApi_backend()
    JuMag.backend[] = oneAPI.oneAPIBackend()
    @info("Intel oneAPI backend is used!")
    return nothing
end

function __init__()
    set_oneApi_backend()
end

end
