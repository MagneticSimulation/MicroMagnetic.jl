module SpinDynamics
using Printf

export create_mesh, create_sim,
       init_m0, add_zeeman, add_dmi,
       add_exch, add_anis, add_demag,
       run_until, relax, save_vtk

include("head.jl")

include("mesh.jl")
include("llg.jl")

include("helper.jl")
include("ode.jl")

include("fileio.jl")
include("sim.jl")
include("demag.jl")

include("vtk.jl")

try
	using CUDAdrv
	include("demag_gpu.jl")
catch
    @warn "CUDA is not available!"
end
end
