# O(1)-storage uniform array. A uniformly valued parameter is "a legal array
# representation", not a separate type: kernels keep a single version that reads
# x[I], and an isbits Fill reaches the kernel by value, so f[I] folds into the
# same immediate as a scalar while saving one memory load.
# Deliberately vendored instead of depending on FillArrays.jl — we need exact
# control over kernel-argument behaviour; rename to `UniformArray` if FillArrays
# is ever adopted.
struct Fill{T,N} <: AbstractArray{T,N}
    value::T
    size::NTuple{N,Int}
end

Fill(v::T, n::Integer) where {T} = Fill{T,1}(v, (Int(n),))

Base.size(f::Fill) = f.size

# Pure computation: no memory access, bounds-free by construction, and the
# compiler folds f[I] into an immediate for isbits element types.
Base.getindex(f::Fill, i::Int) = f.value

Base.IndexStyle(::Type{<:Fill}) = Base.IndexLinear()

Base.sum(f::Fill) = f.value * prod(f.size)

# IO / materialisation path: one O(n) fill instead of per-element getindex.
Base.copyto!(dest::Array, f::Fill) = (fill!(dest, f.value); dest)

Base.:(==)(a::Fill, b::Fill) = (a.value == b.value && a.size == b.size)

# Kernels read α (or similar parameters) per spin; a scalar α (e.g. InertialLLG,
# which keeps `alpha::T`) is wrapped in an O(1) Fill so one kernel serves both.
as_param_array(v::AbstractArray, n::Int) = v
as_param_array(v::Number, n::Int) = Fill(v, n)

# Multiply a per-spin parameter by a scalar, preserving O(1) Fill storage and the
# element type (a plain `v * s` would promote a Fill{Float32} to Fill{Float64}).
scale_param(v::Fill, s) = Fill(convert(eltype(v), v.value * s), length(v))
scale_param(v::AbstractArray, s) = v .* s
