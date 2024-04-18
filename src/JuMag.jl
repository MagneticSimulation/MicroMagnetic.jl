#__precompile__()
module JuMag

using Printf
using KernelAbstractions

const Float = Ref(Float64)
"""
    set_float(x=Float64)
Set the floating-point type for JuMag. 
Defaults to Float64. Use Float32 if single-precision computation is required.
"""
function set_float(x::Type{<:AbstractFloat}=Float64)
    Float[] = x
end

const groupsize = Ref(512)
function set_groupsize(x=n)
    groupsize[] = x
end

const default_backend = Backend[CPU()]
const all_backends = Backend[CPU(), CPU(), CPU(), CPU()]

export set_backend
export set_backend_
@doc raw"""
    set_backend(backend="cuda")
Set the backend of JuMag. Options, hardwares and the corresponding backends are shown as follows: 

| Option                 | Hardware            | Backend                  |
| :---------------------- | :------------------- | :------------------------ |
| "cpu"                  | CPU                 | `KernelAbstractions.CPU()` |
| "cuda" or "nvidia"     | NVIDIA GPU          | `CUDA.CUDABackend()`     |
| "amd"  or "roc"        | AMD GPU             | `AMDGPU.ROCBackend()`    |
| "oneAPI" or "intel"    | Intel GPU           | `oneAPI.oneAPIBackend()` |
| "metal"  or "apple"    | Apple GPU           | `Metal.MetalBackend()`   |

"""
function set_backend(x="cuda")
    backend_names = ["CUDA", "AMDGPU", "oneAPI", "Metal"]
    card_id = 0
    if x == "cuda" || x == "nvidia"
        card_id = 1
    elseif x == "amd" || x == "roc"
        card_id = 2
    elseif x == "oneAPI" || x == "intel"
        card_id = 3
    elseif x == "metal" || x == "apple"
        card_id = 4
    end

    if card_id > 0
        default_backend[] = all_backends[card_id]
        backend_name = backend_names[card_id]
        if Base.find_package(backend_name) === nothing
            @info(@sprintf("Please install %s.jl!", backend_name))
            return
        end

        if default_backend[] == CPU()
            @info(@sprintf("Please import %s!", backend_name))
            return
        end
    else
        default_backend[] = CPU()
    end

    @info(@sprintf("Switch the backend to %s", default_backend[]))
end

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

export using_gpu
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

include("const.jl")
include("micro/mesh.jl")
include("head.jl")
include("util.jl")
include("csg.jl")
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
include("ovf2.jl")
include("vtk.jl")

include("atomistic/mesh.jl")
include("atomistic/head.jl")
include("atomistic/field.jl")
include("atomistic/sim.jl")
include("atomistic/kernels.jl")
include("atomistic/demag.jl")

include("neb/neb.jl")
include("neb/math.jl")
include("neb/neb_kernels.jl")

function __init__() end

end #module
