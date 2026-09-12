# Adaptive mesh refinement (block-structured composite grids) for
# micromagnetics, following García-Cervera & Roma, IEEE Trans. Magn. 42(6),
# 1648 (2006). Public entry point: `AMRSim`.
include("tensor2box.jl")
include("structure.jl")
include("sync.jl")
include("clustering.jl")
include("demag.jl")
include("step.jl")
include("io.jl")
