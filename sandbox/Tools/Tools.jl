module Tools

    using JuMag
    path = dirname(pathof(JuMag))*"/Tools/"
    include(path*"plot.jl")
    include(path*"load_image.jl")
    include(path*"mfm.jl")
    include(path*"xray.jl")
end