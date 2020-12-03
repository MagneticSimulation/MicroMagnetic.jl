#__precompile__()
module JuMC3DTri

using Printf
using CUDA
# using CuArrays
# using CUDA.CUDAdrv

const _cuda_using_double = Ref(false)
const _cuda_available = Ref(true)
# const _mpi_available = Ref(true)


export mu_0, mu_B, k_B, c_e, eV, meV, m_e, g_e, h_bar, gamma, mu_s_1, h_bar_gamma, mT

function cuda_using_double(flag = true)
   _cuda_using_double[] = flag
   return nothing
end

#_cuda_available[] = CuArrays.functional()

if !_cuda_available.x
    @warn "CuArrays does not work!"
end


if _cuda_available.x
    include("mesh.jl")
    include("head.jl")
    include("mc.jl")
    include("util.jl")
    include("mc_kernel.jl")
    include("ovf2.jl")
    include("vtk.jl")
    export MonteCarlo,
           set_exch,
           set_zeeman,
           set_anis_6fold2d,
           init_m0,
           run_step_triangular,
           calcEnAndStd,
           run_sim,
           TriMesh3DGPU,
           save_ovf,read_ovf,
           save_npy,
           save_vtk
end

function __init__()
    _cuda_available[] = CUDA.functional()
    if !_cuda_available.x
        @warn "CUDA is not available!"
    end

end

end #module
