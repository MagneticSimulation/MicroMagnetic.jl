using Documenter
using DemoCards
using JuMag
#using DocumenterTools: Themes

tutorials, tutorials_cb = makedemos("tutorials")

format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    assets=["assets/init.js"]
)

PAGES = ["index.md",
         tutorials,
        "equations.md",
        "notes.md",
        "developer.md",
        "functions.md",
        "questions.md"]

makedocs(
    sitename = "JuMag.jl",
    format = format,
    modules = [JuMag],
    pages = PAGES,
    highlightsig = true
)

tutorials_cb()

deploydocs(
    #deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/ww1g11/JuMag.jl.git"
)
