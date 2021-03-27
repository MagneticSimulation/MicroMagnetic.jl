#__precompile__()
module JuMag

using Printf
using CUDA

const _cuda_using_double = Ref(false)
const _cuda_available = Ref(true)
const _mpi_available = Ref(true)
const _pycall_available = Ref(true)

export init_m0,
       init_m0_random,
       init_m0_skyrmion,init_m0_skyrmion_lattice,
       add_zeeman,
       add_dmi,
       add_exch, add_exch_kagome,
       add_anis, add_cubic_anis,
       add_demag, add_exch_rkky, add_anis_kagome, add_anis_kagome_6fold,
       update_zeeman,update_anis, add_exch_vector,
       run_until, advance_step, relax,
       save_vtk, save_m, ovf2vtk,save_vtk_points,
       FDMesh,
       set_pinning, set_Ms,
       set_Ms_cylindrical,
       Sim,
       CubicMesh,
       TriangularMesh,
       set_mu_s, set_mu_s_kagome,
       set_ux, set_uy, set_uz,
       write_data,
       create_box,create_cylinder,
       compute_system_energy,
       compute_skyrmion_number,
       compute_winding_number_3d,
       winding_number_3d,
       compute_guiding_centre, set_aj,
       NEB,
       interpolate_m,save_ovf,read_ovf,
       fftfreq,
       ex,ey,ez

export mu_0, mu_B, k_B, c_e, eV, meV, m_e, g_e, h_bar, gamma, mu_s_1, h_bar_gamma, mT, Gauss
export ex,ey,ez

function cuda_using_double(flag = true)
   _cuda_using_double[] = flag
   return nothing
end



#this function is copied from CUDA (v.1.7.3) which is gone in the new version
#In the future, we will turn to CUDA since CuArray, CUDAnative are deprecated
function cudims(n::Integer)
  threads = min(n, 256)
  return ceil(Int, n / threads), threads
end

cudims(a::AbstractArray) = cudims(length(a))

include("const.jl")
include("mesh.jl")
include("geometry.jl")
include("head.jl")
include("util.jl")
include("driver.jl")
include("sd.jl")
include("llg/llg.jl")
include("llg/llg_cay.jl")
include("field.jl")
include("helper.jl")
include("integrator/heun.jl")
include("integrator/rk.jl")
include("integrator/dopri5.jl")
include("integrator/dopri5_cayley.jl")
include("fileio.jl")
include("sim.jl")
include("demag.jl")
include("vtk.jl")
include("neb/neb.jl")
include("neb/neb_sd.jl")
include("neb/neb_llg.jl")
include("ovf2.jl")
include("init_m.jl")

#_cuda_available[] = CUDA.functional()

if !_cuda_available.x
    @warn "CUDA does not work!"
end


if _cuda_available.x
    include("cuda/head.jl")
    include("cuda/mesh.jl")
    include("atomistic/head.jl")
    include("atomistic/mesh.jl")
    include("atomistic/sim.jl")
    include("atomistic/field.jl")
    include("cuda/driver.jl")
    include("cuda/sim.jl")
    include("llg/llg_cuda.jl")
    include("llg/llg_cay_cuda.jl")
    include("cuda/util.jl")
    include("cuda/kernels.jl")
    include("cuda/field.jl")
    include("cuda/demag_kernel.jl")
    include("cuda/demag.jl")
    include("cuda/sd.jl")
    include("cuda/vtk.jl")
    include("cuda/ovf2.jl")
    include("neb/neb_cuda.jl")
    include("neb/neb_kernels.jl")
    include("integrator/heun_cuda.jl")
    include("integrator/rk_cuda.jl")
    include("integrator/dopri5_cuda.jl")
    include("integrator/dopri5_cayley_cuda.jl")
    include("mc/mc_helper.jl")
    include("mc/mc.jl")
    include("mc/mc_kernel.jl")
    export FDMeshGPU,
           CubicMeshGPU,
           TriangularMeshGPU,
           MonteCarlo,
           add_thermal_noise,
           set_shape,
           set_shape_to_kagome,
           run_step,
           run_sim,add_exch_anis,
           add_demag_gpu,
           NEB_GPU
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
        CUDA.device!(MPI.Comm_rank(comm) % length(devices()))
    end
end

try
  using PyCall
catch
    _pycall_available[] = false
end

if _pycall_available.x
  include("tools/projection.jl")
  #include("tools/plot.jl")
  include("tools/ltem.jl")
  include("tools/mfm.jl")
  include("tools/xray.jl")
  export vector_field_projection, get_magnetic_phase, OVF2LTEM,OVF2MFM,OVF2XRAY, 
          plotOVF, plot_ovf_projection
end

function __init__()
    _cuda_available[] = CUDA.functional()
    if !_cuda_available.x
        @warn "CUDA is not available!"
    end
    if !_mpi_available.x
        @warn "MPI is not available!"
    end
    if !_pycall_available.x
        @warn "PyCall is not available!"
    end
end

end #module
