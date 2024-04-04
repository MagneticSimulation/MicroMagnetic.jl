#__precompile__()
module JuMag

using Printf
using KernelAbstractions

const single_precision = Ref(false)
"""
    using_single_precision(false)

Set the floating-point precision for simulation, 'true' for single and 'false' for double precision.
By default the double precision will be used. 
"""
function using_single_precision(flag=false)
    single_precision[] = flag
    return nothing
end

const backend = Backend[CPU()]

include("const.jl")
include("micro/mesh.jl")
include("micro/geometry.jl")
include("micro/head.jl")
include("util.jl")
include("micro/driver.jl")
include("micro/sd.jl")
include("llg/llg.jl")
include("micro/kernels.jl")
include("micro/field.jl")
include("micro/add_field.jl")
include("helper.jl")
include("integrator/heun.jl")
include("integrator/rk.jl")
include("integrator/dopri5.jl")
include("fileio.jl")
include("micro/sim.jl")
include("init_m.jl")

end #module
