using Documenter
using JuMag

makedocs(
    sitename = "JuMag.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [JuMag],
    pages = Any[
        "index.md",
        "tutorial.md",
        "notes.md",
        "functions.md"
        ]
)

deploydocs(
    repo = "github.com/ww1g11/JuMag.jl"
)
