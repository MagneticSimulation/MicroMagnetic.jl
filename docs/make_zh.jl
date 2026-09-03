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

# Only translated pages are listed; `pagesonly = true` keeps the rest out of the
# build (and avoids re-running simulations on untranslated example pages).
PAGES = ["首页" => "index.md",
         "手册" => ["install.md", "gui.md", "docker.md", "basics.md", "units.md", "fem.md", "equations.md", "contrib.md"],
         "API" => ["api.md", "api_dev.md"]]

# DOCS_DRAFT=true skips executing example blocks (layout-only iteration runs).
draft = get(ENV, "DOCS_DRAFT", "false") == "true"
# Exported for fillMissingMedia() in src_zh/.vitepress/config.mts (TODO.md #3, Option B).
ENV["DOCUMENTER_MD_ROOT"] = abspath(joinpath(@__DIR__, "build_zh", ".documenter"))

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
    source = "src_zh",
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
