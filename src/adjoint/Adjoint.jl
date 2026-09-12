# Adjoint (VJP) layer.
#
# Hand-written Float64 VJP kernels + linearization for inverse design.  This
# file is pure glue: no module wrapper (the kernels/launchers live in the
# MicroMagnetic namespace next to their production counterparts), and zero
# edits to src/micro, src/eigen, src/atomistic — the field-operator VJP reuses
# the existing `effective_field` launchers and the adjoint of LLGJacOperator is
# added as new mul! methods for its lazy Adjoint/Transpose wrappers.
include("kernels.jl")
include("linearize.jl")

# Steady-state adjoint solver: solves (Df)ᵀλ = P_t∇g at a frozen
# steady state via :pseudotime / :krylov backends, against the VJP callback
# interface (vjp!) — decoupled from the kernels above; apply_KT! from
# linearize.jl plugs into solve_steady_adjoint's vjp! argument.
include("steady.jl")

# Scenario objects: each wires the layer above into the
# run_forward!/gradient!/set_design! user protocol.
include("freqmatch.jl")    # eigenvalue adjoint, FMR matching
