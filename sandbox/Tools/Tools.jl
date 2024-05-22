module Tools

    using MicroMagnetic
    path = dirname(pathof(MicroMagnetic))*"/Tools/"
    include(path*"plot.jl")
    include(path*"load_image.jl")
    include(path*"mfm.jl")
    include(path*"xray.jl")
end