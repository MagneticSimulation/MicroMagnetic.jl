__precompile__()

module JuMag
using Printf

export create_mesh, create_sim,
       init_m0, add_zeeman, add_dmi,
       add_exch, add_anis, add_demag,
	   run_until, relax,
	   save_vtk, FDMesh, set_Ms, Sim,
	   CubicMesh, set_mu_s,
	   set_ux

const _cuda_using_double = Ref(false)
const _cuda_available = Ref(true)

function cuda_using_double(flag = true)
   _cuda_using_double[] = flag
   return nothing
end

include("head.jl")
include("util.jl")

include("mesh.jl")
include("driver.jl")
include("llg.jl")
include("runge_kutta.jl")

include("helper.jl")
include("ode.jl")

include("fileio.jl")
include("sim.jl")
include("demag.jl")
include("vtk.jl")

try
	using CuArrays, CUDAnative #CUDAdrv
    #CuArrays.allowscalar(false) TODO: where should it be?
catch
    _cuda_available[] = false
    @warn "CUDA is not available!"
end

if _cuda_available.x
    include("cuda/head.jl")
    include("cuda/driver.jl")
    include("cuda/sim.jl")
    include("cuda/ode.jl")
    include("cuda/llg.jl")
    include("cuda/util.jl")
    include("cuda/field.jl")
	include("cuda/demag.jl")
    export FDMeshGPU
end


end
