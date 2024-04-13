abstract type Mesh end

export FDMesh

struct FDMesh{T} <: Mesh
    dx::T
    dy::T
    dz::T
    nx::Int64
    ny::Int64
    nz::Int64
    xperiodic::Bool
    yperiodic::Bool
    zperiodic::Bool
    n_total::Int64
    volume::T
    ngbs::AbstractArray{Int32, 2}
end


@inline function index(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64)
    if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
        return -1
    end
    return (k - 1) * nx * ny + (j - 1) * nx + i
end

function indexpbc(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64,
    xperiodic::Bool, yperiodic::Bool, zperiodic::Bool)

    xperiodic && (i < 1) && (i += nx)
    xperiodic && (i > nx) && (i -= nx)

    yperiodic && (j < 1) && (j += ny)
    yperiodic && (j > ny) && (j -= ny)

    zperiodic && (k < 1) && (k += nz)
    zperiodic && (k > nz) && (k -= nz)

    return index(i, j, k, nx, ny, nz)
end


"""
    FDMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")

Create a FDMesh for given parameters. `pbc` could be any combination of "x", "y" and "z".
"""
function FDMesh(; dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")

    ngbs = zeros(Int32, 6, nx * ny * nz)

    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false
    zperiodic = 'z' in pbc ? true : false

    for k = 1:nz, j = 1:ny, i = 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i - 1, j, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
        ngbs[2, id] = indexpbc(i + 1, j, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
        ngbs[3, id] = indexpbc(i, j - 1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
        ngbs[4, id] = indexpbc(i, j + 1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
        ngbs[5, id] = indexpbc(i, j, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
        ngbs[6, id] = indexpbc(i, j, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    end
    volume = dx * dy * dz
    n_total = nx * ny * nz
    if default_backend[] == CPU()
        return FDMesh(dx, dy, dz, nx, ny, nz, xperiodic, yperiodic, zperiodic, n_total, volume, ngbs)
    end

    ngbs_kb = kernel_array(ngbs)
    return FDMesh(dx, dy, dz, nx, ny, nz, xperiodic, yperiodic, zperiodic, n_total, volume, ngbs_kb)
end


@inline function _x_plus_one(i::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, xperiodic::Bool)::Int64
    if i < nx || xperiodic
        return (i == nx) ? index + 1 - nx : index + 1
    end
    return -1
end

@inline function _x_minus_one(i::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, xperiodic::Bool)::Int64
    if i > 1 || xperiodic
        return (i == 1) ? index - 1 + nx : index - 1
    end
    return -1
end

@inline function _y_plus_one(j::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, yperiodic::Bool)::Int64
    if j < ny || yperiodic
        return (j == ny) ? index + nx - nx * ny : index + nx
    end
    return -1
end

@inline function _y_minus_one(j::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, yperiodic::Bool)::Int64
    if j > 1 || yperiodic
        return (j == 1) ? index - nx + nx * ny : index - nx
    end
    return -1
end

@inline function _z_plus_one(k::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, zperiodic::Bool)::Int64
    if k < nz || zperiodic
        return (k == nz) ? index + nx * ny * (1 - nz) : index + nx * ny
    end
    return -1
end

@inline function _z_minus_one(k::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, zperiodic::Bool)::Int64
    if k > 1 || zperiodic
        return (k == 1) ? index + nx * ny * (nz - 1) : index - nx * ny
    end
    return -1
end
