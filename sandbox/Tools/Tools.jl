module Tools

    using NuMag
    path = dirname(pathof(NuMag))*"/Tools/"
    include(path*"plot.jl")
    include(path*"load_image.jl")
    include(path*"mfm.jl")
    include(path*"xray.jl")
end