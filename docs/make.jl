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

deploydocs(
    repo = "github.com/JuliaParallel/MPI.jl.git"
)
