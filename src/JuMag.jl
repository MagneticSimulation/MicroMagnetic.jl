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

end #module
