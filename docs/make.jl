using Documenter
using DemoCards
using MicroMagnetic
#using DocumenterTools: Themes

tutorials, tutorials_cb = makedemos("tutorials")

format = Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true",
                         size_threshold=83886080, assets=["assets/init.js"])

PAGES = ["Home" => "index.md", "basics.md",  "units.md", tutorials, "equations.md", "api.md",
         "api_dev.md", "contrib.md"]

if "warnonly=true" in ARGS
    @info("Option warnonly=true is enabled.")
    makedocs(; sitename="MicroMagnetic.jl", format=format, modules=[MicroMagnetic],
             pages=PAGES, highlightsig=true, checkdocs=:none, warnonly=true)
else
    makedocs(; sitename="MicroMagnetic.jl", format=format, modules=[MicroMagnetic],
             pages=PAGES, highlightsig=true, checkdocs=:none)
end

tutorials_cb()

deploydocs(; repo="github.com/MagneticSimulation/MicroMagnetic.jl")
