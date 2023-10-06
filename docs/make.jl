using Documenter
using DemoCards
using JuMag
#using DocumenterTools: Themes

tutorials, tutorials_cb = makedemos("tutorials")
examples, examples_cb = makedemos("examples")

format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    size_threshold = 83886080,
    assets=["assets/init.js"]
)

PAGES = ["index.md",
        tutorials,
        "equations.md",
         examples,
        "functions.md",
        "notes.md",
        "developer.md",
        "questions.md"]

if "warnonly=true" in ARGS
    @info("Option warnonly=true is enabled.")
    makedocs(
        sitename = "JuMag.jl",
        format = format,
        modules = [JuMag],
        pages = PAGES,
        highlightsig = true,
        checkdocs = :none,
        warnonly=true
    )
else
    makedocs(
        sitename = "JuMag.jl",
        format = format,
        modules = [JuMag],
        pages = PAGES,
        highlightsig = true,
        checkdocs = :none
    )
end

tutorials_cb()
examples_cb()

deploydocs(
    #deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/ww1g11/JuMag.jl.git"
)

