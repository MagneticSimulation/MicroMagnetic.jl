using Documenter
using JuMag

makedocs(
    sitename = "JuMag.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [JuMag],
    pages = Any[
        "notes.md",
        ]
)

#deploydocs(
#    repo = "https://github.com/ww1g11/JuMag.jl"
#)
