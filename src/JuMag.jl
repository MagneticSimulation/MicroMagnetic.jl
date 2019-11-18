#__precompile__()
module JuMag

using Printf
using CUDAnative
using CuArrays
using CuArrays.CUDAdrv

const _cuda_using_double = Ref(false)
const _cuda_available = Ref(true)
const _mpi_available = Ref(true)

export init_m0,
       add_zeeman,
       add_dmi,
       add_exch, add_anis, add_cubic_anis,
       add_demag, add_exch_rkky,
       update_zeeman,update_anis,add_exch_vector,
       run_until, advance_step, relax,
       save_vtk, save_m, ovf2vtk,
       FDMesh, set_Ms, Sim,
       CubicMesh, set_mu_s,
       set_ux, write_data,
       compute_system_energy,
       compute_skyrmion_number,
       compute_winding_number_3d,
       compute_guiding_centre, set_aj,
       NEB,
       interpolate_m,save_ovf,read_ovf

export mu_0, mu_B, k_B, c_e, eV, meV, m_e, g_e, h_bar, gamma, mu_s_1, h_bar_gamma, mT

function cuda_using_double(flag = true)
   _cuda_using_double[] = flag
   return nothing
end

include("const.jl")
include("mesh.jl")
include("head.jl")
include("util.jl")
include("driver.jl")
include("sd.jl")
include("llg.jl")
include("rk.jl")
include("dopri5.jl")
include("helper.jl")
include("ode.jl")
include("fileio.jl")
include("sim.jl")
include("demag.jl")
include("vtk.jl")
include("neb.jl")
include("neb_sd.jl")
include("neb_llg.jl")
include("ovf2.jl")

#_cuda_available[] = CuArrays.functional()

if !_cuda_available.x
    @warn "CuArrays does not work!"
end


if _cuda_available.x
    include("cuda/head.jl")
    include("cuda/mesh.jl")
    include("cuda/driver.jl")
    include("cuda/sim.jl")
    include("cuda/ode.jl")
    include("cuda/llg.jl")
    include("cuda/util.jl")
    include("cuda/kernels.jl")
    include("cuda/field.jl")
    include("cuda/demag_kernel.jl")
    include("cuda/demag.jl")
    include("cuda/sd.jl")
    include("cuda/mc.jl")
    include("cuda/mc_kernel.jl")
    include("cuda/vtk.jl")
    include("cuda/ovf2.jl")
    include("mc/mc.jl")
    include("mc/mc_kernel.jl")
    export FDMeshGPU,
           CubicMeshGPU,
           TriangularMesh,
           MonteCarlo,
           MonteCarloNew,
           run_sim,
           add_demag_gpu
end

try
	using MPI
catch
    _mpi_available[] = false
end

if _mpi_available.x && _cuda_available.x
    include("cuda_mpi/neb.jl")
    include("cuda_mpi/dopri5.jl")
    include("cuda_mpi/neb_kernels.jl")
    export NEB_MPI

    function using_multiple_gpus()
        if !MPI.Initialized()
            MPI.Init()
        end
        comm = MPI.COMM_WORLD
        CUDAnative.device!(MPI.Comm_rank(comm) % length(devices()))
    end
end

function __init__()
    _cuda_available[] = CuArrays.functional()
    if !_cuda_available.x
        @warn "CUDA is not available!"
    end
    if !_mpi_available.x
        @warn "MPI is not available!"
    end
end

end #module
