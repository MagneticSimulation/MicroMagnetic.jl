using Documenter
using MicroMagnetic
using DocumenterVitepress

Atomistic = ["Magnetic skyrmion" => "atomistic/skyrmion.md",
             "Skyrmion lattice" => "atomistic/skyrmion_lattice.md",
             "AFM skyrmion" => "atomistic/skyrmion_afm.md",
             "Stochastic LLG" => "atomistic/sllg.md",
             "Phase diagram" => "atomistic/phase_diagram.md"];

Micromagnetic = ["Nanobar" => "micromagnetics/nanobar.md",
                 "Magnetic vortex" => "micromagnetics/vortex.md",
                 "Standard Problem 4 (sim_with)" => "micromagnetics/std4_sim_with.md",
                 "Standard Problem 4" => "micromagnetics/std4.md",
                 "Standard Problem 5 (sim_with)" => "micromagnetics/std5_sim_with.md",
                 "Standard Problem 5" => "micromagnetics/std5.md",
                 "Skyrmion dynamics STT" => "micromagnetics/skyrmion_stt.md",
                 "Stoner–Wohlfarth model" => "micromagnetics/stoner_wohlfarth.md",
                 "Dynamical susceptibility" => "micromagnetics/chi.md"]

Miscellaneous = ["Skyrmion Phase (Monte Carlo)" => "monte_carlo/skyrmion.md",
                 "M-T curve (Monte Carlo)" => "monte_carlo/M_T_curve.md",
                 "Skyrmion collapse (NEB)" => "neb/neb_skx.md"]

PAGES = ["Home" => "index.md",
         "Manual" => ["install.md", "basics.md", "units.md", "equations.md"],
         "Atomistic" => Atomistic, 
         "Micromagnetic" => Micromagnetic, 
         "Miscellaneous" => Miscellaneous, 
         "api.md", "api_dev.md",
         "contrib.md"]


makedocs(; 
    sitename = "MicroMagnetic.jl", 
    modules = [MicroMagnetic],
    warnonly = true,
    checkdocs=:all,
    format= MarkdownVitepress(; repo="github.com/MagneticSimulation/MicroMagnetic.jl",
                           devbranch="master", devurl="dev"),
    draft = false,
    source = "src",
    build = "build",
    pages = PAGES
)

#deploydocs(; repo="github.com/MagneticSimulation/MicroMagnetic.jl")
DocumenterVitepress.deploydocs(; repo="github.com/MagneticSimulation/MicroMagnetic.jl",
                               target="build")
