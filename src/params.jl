export set_region, region_map

"""
    make_param(T, v, mesh, n_total; scale=1)

Materialise the parameter `v` (a `NumberOrArrayOrFunction`) into an n_total-length
read-only array of element type `T`. This is the single discrimination point for
"spatial or not": a uniform number becomes an O(1) `Fill`, a
function or array is materialised into a dense backend array through the existing
`init_scalar!` path. `scale` (e.g. `mu_0` for `set_Ms`) is applied exactly as the
old "store, then `. *= scale`" two-step did, keeping results bit-identical.

For non-isbits element types (the symbolic `AbstractFloat` mode) the number falls
back to dense materialisation, because `Fill{AbstractFloat}` is not isbits and
cannot reach GPU kernels.
"""
function make_param(::Type{T}, v::NumberOrArrayOrFunction, mesh::Mesh, n_total::Int;
                    scale=1) where T
    if v isa Number && isbitstype(T)
        # round to T first, then round again after scaling: exactly what the dense
        # "store, then .*= scale" path computes. Skipping the final
        # T() would promote to Float64 for scale=mu_0 and silently materialise a
        # dense CPU array at the field assignment (eltype mismatch).
        v = T(v)
        v = scale == 1 ? v : T(v * scale)
        return Fill(v, n_total)
    end
    a_host = zeros(T, n_total)
    init_scalar!(a_host, mesh, v)
    if scale != 1
        a_host .*= scale
    end
    if default_backend[] == CPU()
        return a_host
    end
    # a CPU array cannot be broadcast into a GPU array, so materialise on the host
    # and make one copy to the backend (the old zeros → init_scalar! → copy pattern)
    a = create_zeros(T, n_total)
    copyto!(a, a_host)
    return a
end

"""
    make_vector_param(T, v, mesh, n_total) -> (Hx, Hy, Hz)

Three-component variant of [`make_param`](@ref): a 3-tuple input maps component-wise
(numbers become O(1) `Fill`s, functions dense arrays), while an n_total-array or a
function returning 3-tuples is materialised into three dense component segments.
"""
function make_vector_param(::Type{T}, v::TupleOrArrayOrFunction, mesh::Mesh,
                           n_total::Int) where T
    if v isa Tuple && length(v) == 3
        return (make_param(T, v[1], mesh, n_total),
                make_param(T, v[2], mesh, n_total),
                make_param(T, v[3], mesh, n_total))
    end
    f = zeros(T, 3 * n_total)
    init_vector!(f, mesh, v)
    # the 3N layout is spin-interleaved (mx1,my1,mz1,mx2,...): split per component
    # on the host with strided copies before moving to the backend
    Hx = kernel_array(f[1:3:3 * n_total])
    Hy = kernel_array(f[2:3:3 * n_total])
    Hz = kernel_array(f[3:3:3 * n_total])
    return (Hx, Hy, Hz)
end

"""
    set_region(mesh::FDMesh, region_id::Int, shape::CSGShape)
    set_region(mesh::FDMesh, shape::CSGShape, region_id=1)

Set the region ID for cells inside the specified shape.

Parameters:
- mesh: The mesh to set regions for
- region_id: The ID to assign to cells inside the shape (default: 1)
- shape: The shape defining the region

# Examples
```julia
# Create a finite difference mesh
mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=10, ny=10, nz=10)

# Create a spherical shape centered at (10e-9, 10e-9, 0) with radius 10e-9
sphere = Sphere(radius=10e-9, center=(10e-9, 10e-9, 0))

# Option 1: With positional region_id (backward compatible)
set_region(mesh, 1, sphere)
# or
set_region(mesh, sphere, 1)

# Option 2: With keyword region_id (more readable)
set_region(mesh, sphere; region_id=1)
```
"""
function set_region(mesh::FDMesh, region_id::Int, shape::CSGShape)
    a = isa(mesh.regions, Array) ? mesh.regions : Array(mesh.regions)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        x = mesh.x0 + (i - 0.5)*dx
        y = mesh.y0 + (j - 0.5)*dy
        z = mesh.z0 + (k - 0.5)*dz
        
        if inside(shape, (x, y, z))
            a[id] = region_id
        end
    end
    
    isa(mesh.regions, Array) || copyto!(mesh.regions, a)
    mesh.layout_version[] += 1
    return true
end

# Overloaded method with shape first and optional region_id
function set_region(mesh::FDMesh, shape::CSGShape; region_id::Int=1)
    return set_region(mesh, region_id, shape)
end
# This allows calling as set_region(mesh, sphere, 1) or set_region(mesh, sphere, id=1)
function set_region(mesh::FDMesh, shape::CSGShape, region_id::Int)
    return set_region(mesh, region_id, shape)
end

"""
    region_map(mapping...; default=0.0)

Create a region mapping function that maps region IDs to values

Parameters
- `mapping::Pair{Int,<:Number}`: Mapping pairs from region_id to value
- `default::Number`: Default value (when region_id is not in the mapping)

Return Value
A function `f(region_id)` that returns the corresponding value according to the mapping

Examples
```julia
# Simple mapping for two regions
Ms_func = region_map(1 => 8.6e5)  # region 1: 8.6e5, others: 0.0

# Mapping for multiple regions
exch_func = region_map(
    1 => 1.3e-11,
    2 => 0.8e-11,
    3 => 0.5e-11
)

# Direct usage
set_Ms(sim, region_map(1 => 8.6e5))
```
"""
function region_map(mapping::Pair{Int,<:Number}...; default=0.0)
    dict = Dict{Int,Float64}()
    for (region_id, value) in mapping
        dict[region_id] = Float64(value)
    end
    return function (region_id::Integer)
        return get(dict, region_id, default)
    end
end
