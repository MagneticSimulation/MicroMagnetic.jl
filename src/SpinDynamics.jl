__precompile__()

module SpinDynamics
using Printf

export create_mesh, create_sim,
       init_m0, add_zeeman, add_dmi,
       add_exch, add_anis, add_demag,
	   run_until, relax,
	   save_vtk, FDMesh, set_Ms, Sim,
	   CubicMesh, set_mu_s

const FloatGPU = Float64

include("head.jl")

include("mesh.jl")
include("driver.jl")
include("llg.jl")

include("helper.jl")
include("ode.jl")

include("fileio.jl")
include("sim.jl")
include("demag.jl")

include("vtk.jl")

cuda_available = true
try
	using CuArrays, CUDAnative, CUDAdrv
catch
    cuda_available = false
    @warn "CUDA is not available!"
end

if cuda_available
    include("demag_gpu.jl")
    include("cuda/head.jl")
    include("cuda/driver.jl")
    include("cuda/sim.jl")
    include("cuda/ode.jl")
    include("cuda/util.jl")
    include("cuda/field.jl")
    export add_demag_gpu, FDMeshGPU
end



end
