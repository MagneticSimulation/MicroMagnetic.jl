#GPU smoke entry for the adjoint track.
#
#Standalone run:  julia --project=. test/adjoint/cuda/runtests.jl
#
#T1's adjoint kernel tests slot in here as extra `include`s, so every kernel gets
#verified on CUDA the moment it lands (CPU-first in test/adjoint/, GPU 紧随).
using Test
using MicroMagnetic

#CUDA is a hard requirement for this directory: fail loudly instead of silently
#skipping (the soft-skip convention lives in test/test_utils.jl's test_functions).
Base.find_package("CUDA") === nothing &&
    error("CUDA.jl is not loadable; test/adjoint/cuda requires it (weakdep: `using CUDA` in an environment that has it)")
using CUDA
CUDA.allowscalar(false)

@testset "adjoint/cuda" begin
    include("test_smoke.jl")
end
