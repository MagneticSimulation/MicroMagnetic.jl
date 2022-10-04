using Documenter
using DemoCards
using JuMag
#using DocumenterTools: Themes

examples, examples_cb = makedemos("examples")

format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    assets=["assets/init.js"]
)

PAGES = ["index.md",
    "equations.md",
    examples,
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

examples_cb()

deploydocs(
    #deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/ww1g11/JuMag.jl.git"
)
