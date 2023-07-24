#__precompile__()
module JuMag

using LinearAlgebra
using Printf
using CUDA

function dev_test()
    return "This is a newly added function!"
end

LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

const _cuda_using_double = Ref(false)
const _cuda_available = Ref(true)
const _mpi_available = Ref(true)

export init_m0,
       init_m0_random,
       init_m0_skyrmion, 
       init_m0_skyrmion_lattice,
       add_zeeman,
       add_dmi,
       add_exch, add_exch_kagome, add_next_exch, add_next_next_exch, add_next_next_next_exch,
       add_anis, add_cubic_anis,
       add_demag, add_exch_rkky, add_anis_kagome, add_anis_kagome_6fold,
       update_zeeman,update_anis, add_exch_vector,
       run_until, advance_step, relax,
       save_vtk, save_m, ovf2vtk,save_vtk_points, mag2ovf,
       FDMesh,
       set_pinning, set_Ms,
       set_Ms_cylindrical,
       set_driver,
       Sim,
       create_sim,
       run_sim, 
       SaverItem, DataSaver,init_saver,
       CubicMesh,
       TriangularMesh,
       set_mu_s, set_mu_s_kagome,
       set_ux, set_uy, set_uz,
       write_data,save_sim_data,
       create_box,create_cylinder,
       compute_system_energy,
       compute_skyrmion_number,
       compute_tensor_B,
       compute_winding_number_3d,
       winding_number_3d,
       compute_guiding_centre, set_aj,
       compute_guiding_center,
       effective_field,
       NEB,
       interpolate_m,save_ovf,read_ovf,
       fftfreq,
       ex,ey,ez,
       jdl2png,
       jdl2movie

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
include("add_field.jl")
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
include("movie.jl")

#include("main.jl")


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
    include("atomistic/kernels.jl")
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
    include("atomistic/demag_kernel.jl")
    include("atomistic/demag.jl")
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
           CubicMesh4NGPU,
           SquareMeshGPU,
           TriangularMeshGPU,
           CylindricalTubeMeshGPU,
           MonteCarlo,
           add_thermal_noise,
           set_shape,
           set_shape_to_kagome,
           run_step,
           run_sim,add_exch_anis,
           add_demag_gpu,
           NEB_GPU,
           add_magnetoelectric_laser,
           add_anis_tube
end


include("ltem/interpolation.jl")
include("ltem/tilt.jl")
include("ltem/projection.jl")
include("ltem/ltem.jl")

export compute_magnetic_phase, warp, radon


const _tools_loaded = Ref(false)
function load_tools()
    if !_tools_loaded[]
        module_path = dirname(pathof(JuMag))
        include(module_path * "/Tools/Tools.jl")
        _tools_loaded[] = true
    end
end


function __init__()
    _cuda_available[] = CUDA.functional()
    if !_cuda_available.x
        @warn "CUDA is not available!"
    end
end

end #module
