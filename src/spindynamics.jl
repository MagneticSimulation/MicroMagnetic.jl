module SpinDynamics
using Printf

export create_mesh, create_sim,
       init_m0, add_zeeman,
       add_exch, add_anis, add_demag,
       run_until, relax

include("head.jl")

include("mesh.jl")
include("llg.jl")

include("helper.jl")
include("ode.jl")

include("fileio.jl")
include("sim.jl")
include("demag.jl")

end
