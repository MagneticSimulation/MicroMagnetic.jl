using Documenter
using DemoCards
using MicroMag
#using DocumenterTools: Themes

#tutorials, tutorials_cb = makedemos("tutorials")
#examples, examples_cb = makedemos("examples")

format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    size_threshold = 83886080,
    assets=["assets/init.js"]
)

PAGES = ["Home" =>"index.md",
        "basics.md",
        #tutorials,
        #"equations.md",
        #examples,
        "functions.md",
        "developer.md"]

if "warnonly=true" in ARGS
    @info("Option warnonly=true is enabled.")
    makedocs(
        sitename = "MicroMag.jl",
        format = format,
        modules = [MicroMag],
        pages = PAGES,
        highlightsig = true,
        checkdocs = :none,
        warnonly=true
    )
else
    makedocs(
        sitename = "MicroMag.jl",
        format = format,
        modules = [MicroMag],
        pages = PAGES,
        highlightsig = true,
        checkdocs = :none
    )
end


#tutorials_cb()
#examples_cb()

deploydocs(
    repo = "github.com/ww1g11/MicroMag.jl.git"
)

