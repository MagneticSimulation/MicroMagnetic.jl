using Documenter
using JuMag

makedocs(
    sitename = "JuMag",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [JuMag],
    pages = Any[
        "index.md",
        "notes.md"
        ]
)

deploydocs(
    repo = "github.com/ww1g11/JuMag.jl"
)
