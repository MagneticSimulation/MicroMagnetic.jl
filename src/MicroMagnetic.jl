#__precompile__()
module MicroMagnetic

using Printf
using KernelAbstractions

using TimerOutputs
const timer = TimerOutput()

function reset_timer()
    return reset_timer!(timer)
end
export reset_timer

const Float = Ref(Float64)
"""
    set_precision(x::Type{<:AbstractFloat}=Float64)

Set the precision for MicroMagnetic simulations.

This function allows you to specify the precision of floating-point numbers to be used in MicroMagnetic simulations. 
By default, it sets the precision to `Float64`. If single-precision computation is required, you can specify `Float32`.

# Example
```julia
set_precision(Float32)
```
"""
function set_precision(x::Type{<:AbstractFloat}=Float64)
    return Float[] = x
end
export set_precision

Base.@deprecate set_float(x::Type{<:AbstractFloat}=Float64) set_precision(x)

const groupsize = Ref(512)
function set_groupsize(x)
    return groupsize[] = x
end
export set_groupsize

macro set_backend(platform)
    return esc(quote
                   if Base.find_package($(QuoteNode(platform))) !== nothing
                       using $(Symbol(platform))
                   end
               end)
end

const default_backend = Backend[CPU()]
const all_backends = Backend[CPU(), CPU(), CPU(), CPU()]


export set_backend
@doc raw"""

    set_backend(backend="cuda")

Set the backend of MicroMagnetic. This function allows you to specify the backend for MicroMagnetic simulations. 

The available options and their corresponding hardware and backends are shown below:

| Option                 | Hardware            | Backend                  |
| :--------------------- | :------------------ | :------------------------ |
| "cpu"                  | CPU                 | `KernelAbstractions.CPU()` |
| "cuda" or "nvidia"     | NVIDIA GPU          | `CUDA.CUDABackend()`     |
| "amd" or "roc"         | AMD GPU             | `AMDGPU.ROCBackend()`    |
| "oneAPI" or "intel"    | Intel GPU           | `oneAPI.oneAPIBackend()` |
| "metal" or "apple"     | Apple GPU           | `Metal.MetalBackend()`   |

# Examples

To set the backend to use CUDA (NVIDIA GPU):

```julia
using MicroMagnetic
using CUDA
```

To set the backend to use the CPU/CUDA:
```
set_backend("cpu")
set_backend("cuda")
```

"""
function set_backend(backend="cuda")
    backend_names = ["CUDA", "AMDGPU", "oneAPI", "Metal"]
    card_id = 0
    x = lowercase(backend)
    if x == "cuda" || x == "nvidia"
        card_id = 1
    elseif x == "amd" || x == "roc" || x == "amdgpu"
        card_id = 2
    elseif x == "oneapi" || x == "intel"
        card_id = 3
    elseif x == "metal" || x == "apple"
        card_id = 4
    end

    if card_id > 0
        default_backend[] = all_backends[card_id]
        backend_name = backend_names[card_id]
        if Base.find_package(backend_name) === nothing
            @info(@sprintf("Please install %s.jl!", backend_name))
            return false
        end

        if default_backend[] == CPU()
            @info(@sprintf("Please import %s!", backend_name))
            return false
        end
    else
        default_backend[] = CPU()
    end

    @info(@sprintf("Switch the backend to %s", default_backend[]))
    return true
end


# A global reference to control verbose logging
# Verbose is a Ref that holds a Boolean value indicating whether verbose logging is enabled
const Verbose = Ref(false)

"""
    set_verbose(verbose::Bool=false)

Sets the global verbosity level for logging.

# Arguments
- `verbose::Bool`: If true, enables verbose logging; if false, disables it. Default is false.

# Example
```julia
set_verbose(true)  # Enables verbose logging
set_verbose(false) # Disables verbose logging
```
"""
function set_verbose(verbose::Bool=true) 
    Verbose[] = verbose 
end
export set_verbose

function kernel_array(a::Array)
    A = KernelAbstractions.zeros(default_backend[], eltype(a), size(a))
    copyto!(A, a)
    return A
end

function create_zeros(dims...)
    T = Float[]
    return KernelAbstractions.zeros(default_backend[], T, dims)
end

function create_zeros(::Type{T}, dims...) where {T}
    return KernelAbstractions.zeros(default_backend[], T, dims)
end

function create_ones(::Type{T}, dims...) where {T}
    return KernelAbstractions.ones(default_backend[], T, dims)
end

export @using_gpu
macro using_gpu()
    quote
        try
            using CUDA
        catch
        end
        try
            using AMDGPU
        catch
        end
        try
            using oneAPI
        catch
        end
        try
            using Metal
        catch
        end
    end
end

function ovf2png() end
function plot_m() end
function jld2movie() end
function dynamic_matrix() end

export ovf2png, plot_m, jld2movie, dynamic_matrix

include("const.jl")
include("micro/mesh.jl")
include("head.jl")
include("util.jl")
include("csg.jl")
include("driver.jl")
include("sd.jl")
include("llg/llg.jl")
include("llg/llg_cayley.jl")
include("micro/kernels.jl")
include("micro/field.jl")
include("micro/demag.jl")
include("micro/demag_direct.jl")
include("micro/add_field.jl")
include("helper.jl")
include("integrator/heun.jl")
include("integrator/rk.jl")
include("integrator/rk_cayley.jl")
include("integrator/dopri5.jl")
include("integrator/dopri5_cayley.jl")
include("fileio.jl")
include("sim.jl")
include("init_m.jl")
include("ovf2.jl")
include("vtk.jl")
include("voronoi.jl")

include("atomistic/mesh.jl")
include("atomistic/head.jl")
include("atomistic/field.jl")
include("atomistic/sim.jl")
include("atomistic/kernels.jl")
include("atomistic/demag.jl")

include("mc/kernels.jl")
include("mc/mc_helper.jl")
include("mc/mc.jl")

include("neb/neb.jl")
include("neb/math.jl")
include("neb/neb_kernels.jl")

include("eigen/util.jl")
include("eigen/eigen.jl")

include("tools/ltem.jl")

include("fem/mesh.jl")
include("fem/util.jl")

function __init__() end

end #module
