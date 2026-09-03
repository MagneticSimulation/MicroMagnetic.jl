# docs/TODO.md #7: translation drift check (warn only, never fails the build).
#
# For every translated page docs/src_zh/<path>, compare the last commit that touched
# the English source docs/src/<path>; warn when the source changed after the
# translation was last updated. Run from anywhere:
#
#     julia docs/check_zh_drift.jl
#
using Dates

const REPO = dirname(@__DIR__)
const SRC = joinpath("docs", "src")
const SRC_ZH = joinpath("docs", "src_zh")

function last_commit(rel)
    out = readchomp(`git -C $REPO log -1 --format=%cI -- $(joinpath(SRC_ZH, rel))`)
    zh = isempty(out) ? nothing : DateTime(out[1:19])
    out = readchomp(`git -C $REPO log -1 --format=%cI -- $(joinpath(SRC, rel))`)
    en = isempty(out) ? nothing : DateTime(out[1:19])
    return zh, en
end

function main()
    outdated = String[]
    n_translated = 0
    n_skipped = 0
    for (root, dirs, files) in walkdir(joinpath(REPO, SRC_ZH)), file in files
        endswith(file, ".md") || continue
        rel = relpath(joinpath(root, file), joinpath(REPO, SRC_ZH))
        isfile(joinpath(REPO, SRC, rel)) || continue  # zh-only page, nothing to compare
        n_translated += 1
        zh, en = last_commit(rel)
        if zh === nothing || en === nothing
            n_skipped += 1  # uncommitted file, skip
            continue
        end
        if en > zh
            push!(outdated, rel)
        end
    end

    if isempty(outdated)
        skipped = n_skipped > 0 ? " ($n_skipped uncommitted skipped)" : ""
        println("✅ Translation drift check: all $n_translated translated pages are up to date$skipped.")
    else
        @warn "Translation drift detected: the English source of the following page(s) changed after the last translation update. Consider re-syncing them (this warning does not fail the build)."
        for rel in outdated
            println("  - docs/src_zh/$rel")
        end
    end
end

main()
