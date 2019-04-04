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

try
	using CUDAnative
	include("demag_gpu.jl")
    include("gpu/head.jl")
    export add_demag_gpu, FDMeshGPU
catch
    @warn "CUDA is not available!"
end
end
