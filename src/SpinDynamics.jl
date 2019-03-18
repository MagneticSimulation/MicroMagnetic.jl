module SpinDynamics
using Printf

export create_mesh, create_sim,
       init_m0, add_zeeman, add_dmi,
       add_exch, add_anis, add_demag,
       add_demag_gpu, 
       run_until, relax, save_vtk

include("head.jl")

include("mesh.jl")
include("llg.jl")

include("helper.jl")
include("ode.jl")

include("fileio.jl")
include("sim.jl")
include("demag.jl")
include("demag_gpu.jl")

include("vtk.jl")
end
