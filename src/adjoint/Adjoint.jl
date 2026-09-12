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
