# Chinese documentation build (docs/TODO.md #7).
#
# Deployed to the SAME gh-pages branch as the English site, in a top-level `zh/`
# folder (sibling of `dev/`): https://magneticsimulation.github.io/MicroMagnetic.jl/zh/
# Each deploydocs run only replaces its own devurl folder, so the English `dev/`
# tree and this `zh/` tree never clobber each other. See TODO.md for the spike that
# validated this (and why `dev/zh/` nesting was rejected).
using Documenter
using MicroMagnetic
using DocumenterVitepress
using CairoMakie

# Untranslated pages fall back to their verbatim English source so the Chinese
# site stays complete. Each build assembles the gitignored union dir
# `src_zh_all/`: the English tree (`src/`) as the base, with the translated
# pages (`src_zh/`) overlaid on top — a page is Chinese as soon as
# `src_zh/<page>` exists, nothing else needs to change. Overlaying (instead of
# copying pages one by one) also brings along the per-page data files the
# example blocks read (e.g. `micromagnetics/assets/chi.txt`).
UNION = joinpath(@__DIR__, "src_zh_all")
rm(UNION; recursive = true, force = true)
cp(joinpath(@__DIR__, "src"), UNION)
ZH = joinpath(@__DIR__, "src_zh")
for (root, dirs, files) in walkdir(ZH), f in files
    dst = joinpath(UNION, relpath(joinpath(root, f), ZH))
    mkpath(dirname(dst))
    cp(joinpath(root, f), dst; force = true)
end

PAGES = ["首页" => "index.md",
         "用户手册" => ["install.md", "basics.md", "units.md", "gpu.md", "pbc.md",
                     "dataio.md", "tools.md", "fem.md", "eigen/eigenmodes.md",
                     "gui.md", "docker.md"],
         "物理" => ["equations.md"],
         "原子模型" => ["atomistic/skyrmion.md",
                     "atomistic/saddle_point_search.md",
                     "atomistic/skyrmion_lattice.md",
                     "atomistic/hopfion.md",
                     "atomistic/skyrmion_afm.md",
                     "atomistic/sllg.md",
                     "atomistic/phase_diagram.md"],
         "微磁模型 (FD)" => ["micromagnetics/nanobar.md",
                          "micromagnetics/vortex.md",
                          "micromagnetics/std4.md",
                          "micromagnetics/std4_sim_with.md",
                          "micromagnetics/std5.md",
                          "micromagnetics/std5_sim_with.md",
                          "micromagnetics/stoner_wohlfarth.md",
                          "micromagnetics/chi.md"],
         "微磁模型 (FE)" => ["fem/sphere_demag.md",
                          "fem/hysteresis.md",
                          "fem/rkky.md"],
         "杂项" => ["monte_carlo/skyrmion.md",
                  "monte_carlo/M_T_curve.md",
                  "neb/neb_skx.md"],
         "API" => ["api.md"],
         "开发者" => ["contrib.md"]]

# DOCS_DRAFT=true skips executing example blocks (layout-only iteration runs).
draft = get(ENV, "DOCS_DRAFT", "false") == "true"

# On tag pushes (docs.yml triggers on `tags: '*'`) the English build deploys
# versioned folders (`vX.Y.Z/`, `stable/`). Without pinning, this build's
# deploy_decision would resolve to the same version subfolder and the Chinese
# deploy would overwrite the English versioned docs with the (partial) Chinese
# tree. Pin the subfolder to `zh` so this build only ever owns `zh/`.
zh_deploy_decision = if startswith(get(ENV, "GITHUB_REF", ""), "refs/tags/")
    Documenter.DeployDecision(; all_ok = true,
                              repo = "github.com/MagneticSimulation/MicroMagnetic.jl",
                              branch = "gh-pages", is_preview = false, subfolder = "zh")
else
    nothing
end

makedocs(;
    sitename = "MicroMagnetic.jl",
    modules = [MicroMagnetic,
              isdefined(Base, :get_extension) ? Base.get_extension(MicroMagnetic, :CairoMakieExt) : MicroMagnetic.CairoMakieExt],
    warnonly = true,
    checkdocs=:all,
    format= MarkdownVitepress(; repo="github.com/MagneticSimulation/MicroMagnetic.jl",
                           devbranch="master", devurl="zh", deploy_decision=zh_deploy_decision),
    draft = draft,
    source = "src_zh_all",
    build = "build_zh",
    pagesonly = true,
    pages = PAGES
)

if get(ENV, "CI", "false") == "true"
    # NOTE: do not pass `devurl` here. Documenter.deploydocs must keep its default
    # `devurl="dev"` so that versions.js and the root redirect keep pointing at the
    # English `dev` site after the Chinese deploy rewrites them.
    DocumenterVitepress.deploydocs(; repo="github.com/MagneticSimulation/MicroMagnetic.jl",
                                   target="build_zh")
end
