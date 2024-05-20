module Tools

    using MicroMag
    path = dirname(pathof(MicroMag))*"/Tools/"
    include(path*"plot.jl")
    include(path*"load_image.jl")
    include(path*"mfm.jl")
    include(path*"xray.jl")
end