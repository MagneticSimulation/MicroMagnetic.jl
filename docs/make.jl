using Documenter
using MicroMagnetic
using DocumenterVitepress
using CairoMakie

Atomistic = ["Magnetic skyrmion" => "atomistic/skyrmion.md",
             "Saddle-point search (SPS)" => "atomistic/saddle_point_search.md",
             "Skyrmion lattice" => "atomistic/skyrmion_lattice.md",
             "Magnetic hopfion" => "atomistic/hopfion.md",
             "AFM skyrmion" => "atomistic/skyrmion_afm.md",
             "Stochastic LLG" => "atomistic/sllg.md",
             "Phase diagram" => "atomistic/phase_diagram.md"];

Micromagnetic = ["Nanobar" => "micromagnetics/nanobar.md",
                 "Magnetic vortex" => "micromagnetics/vortex.md",
                 "Standard Problem 4" => "micromagnetics/std4.md",
                 "Standard Problem 4 (sim_with)" => "micromagnetics/std4_sim_with.md",
                 "Standard Problem 5" => "micromagnetics/std5.md",
                 "Standard Problem 5 (sim_with)" => "micromagnetics/std5_sim_with.md",
                 #"Skyrmion dynamics STT" => "micromagnetics/skyrmion_stt.md",
                 "Stoner-Wohlfarth model" => "micromagnetics/stoner_wohlfarth.md",
                 "Dynamical susceptibility" => "micromagnetics/chi.md"]

FE = ["Magnetized Sphere" => "fem/sphere_demag.md",
      "Hysteresis Loop" => "fem/hysteresis.md"]

API = ["api.md", "api_dev.md"]

Miscellaneous = ["Skyrmion Phase (Monte Carlo)" => "monte_carlo/skyrmion.md",
                 "M-T curve (Monte Carlo)" => "monte_carlo/M_T_curve.md",
                 "Skyrmion collapse (NEB)" => "neb/neb_skx.md"]

PAGES = ["Home" => "index.md",
         "Manual" => ["install.md", "gui.md", "docker.md", "basics.md", "units.md", "fem.md", "equations.md", "eigen/eigenmodes.md", "contrib.md"],
         "Atomistic" => Atomistic,
         "Micromagnetics (FD)" => Micromagnetic,
         "Micromagnetics (FE)" => FE,
         "Miscellaneous" => Miscellaneous,
         "API" => API]


# DOCS_DRAFT=true skips executing example blocks (layout-only iteration runs).
draft = get(ENV, "DOCS_DRAFT", "false") == "true"
# Exported for fillMissingMedia() in src/.vitepress/config.mts (docs/TODO.md #3, Option B).
ENV["DOCUMENTER_MD_ROOT"] = abspath(joinpath(@__DIR__, "build", ".documenter"))

makedocs(;
    sitename = "MicroMagnetic.jl",
    modules = [MicroMagnetic,
              isdefined(Base, :get_extension) ? Base.get_extension(MicroMagnetic, :CairoMakieExt) : MicroMagnetic.CairoMakieExt],
    warnonly = true,
    checkdocs=:all,
    format= MarkdownVitepress(; repo="github.com/MagneticSimulation/MicroMagnetic.jl",
                           devbranch="master", devurl="dev"),
    draft = draft,
    source = "src",
    build = "build",
    pagesonly = true,
    pages = PAGES
)

#deploydocs(; repo="github.com/MagneticSimulation/MicroMagnetic.jl")
if get(ENV, "CI", "false") == "true"
    DocumenterVitepress.deploydocs(; repo="github.com/MagneticSimulation/MicroMagnetic.jl",
                                   target="build")
end
