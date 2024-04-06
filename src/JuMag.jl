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
const all_backends = Backend[CPU(), CPU(), CPU(), CPU()]

export set_backend

"""
    set_backend(backend="cuda")
Set the backend of JuMag, options are: 
    - "cpu"                :  CPU
    - "cuda" or "nvidia"   :  NVIDIA GPU
    - "amd"                :  AMD GPU
    - "oneAPI" or "intel"  :  Intel GPU 
    - "metal"  or "apple"  :  Apple GPU
"""
function set_backend(x="cuda")
    backend_names = ["CUDA", "AMDGPU", "oneAPI", "Metal"]
    card_id = 0
    if x == "cuda" || x == "nvidia"
        card_id = 1
    elseif x == "amd" 
        card_id = 2
    elseif x == "oneAPI" ||  x == "intel"
        card_id = 3
    elseif x == "metal"  ||  x == "apple"
        card_id = 4
    end

    if card_id > 0
        backend[] = all_backends[card_id]
        backend_name = backend_names[card_id]
        if Base.find_package(backend_name) === nothing
            @info(@sprintf("Please install %s.jl!", backend_name))
            return 
        end

        if backend[] == CPU()
            @info(@sprintf("Please import %s!", backend_name))
            return
        end
    end

    @info(@sprintf("Switch the backend to %s", backend[]))
    
end

include("const.jl")
include("micro/mesh.jl")
include("micro/geometry.jl")
include("head.jl")
include("util.jl")
include("micro/driver.jl")
include("micro/sd.jl")
include("llg/llg.jl")
include("micro/kernels.jl")
include("micro/field.jl")
include("micro/demag.jl")
include("micro/add_field.jl")
include("helper.jl")
include("integrator/heun.jl")
include("integrator/rk.jl")
include("integrator/dopri5.jl")
include("fileio.jl")
include("sim.jl")
include("init_m.jl")

end #module
