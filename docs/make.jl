using Documenter
#using DemoCards
using JuMag
#using DocumenterTools: Themes

format = Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true")

PAGES = ["index.md",
    "tutorial.md",
    "equations.md",
    "notes.md",
    "functions.md",
    "questions.md"]

makedocs(
    sitename = "JuMag.jl",
    format = format,
    modules = [JuMag],
    pages = PAGES,
    highlightsig = true
)


deploydocs(
    #deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/ww1g11/JuMag.jl.git"
)
