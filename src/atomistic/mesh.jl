abstract type AtomisticMesh <: Mesh end

export CubicMesh, TriangularMesh, CylindricalTubeMesh

struct TriangularMesh <: AtomisticMesh
    dx::Float64
    dy::Float64
    dz::Float64
    nx::Int64
    ny::Int64
    nz::Int64
    n_total::Int64
    n_ngbs::Int64
    n_ngbs2::Int64
    ngbs::AbstractArray{Int32,2}
    ngbs2::AbstractArray{Int32,2}
    xperiodic::Bool
    yperiodic::Bool
    zperiodic::Bool
end

@doc raw"""
    TriangularMesh(; dx=1e-9, nx=1, ny=1, pbc="open")

Create a 2d triangular mesh. The index of the nearest neighbours and the next-nearest 
neighbours are given as follows:

| nearest index | location   | next-nearest index | location |
| :--------: | :------------: | :---------: | :---------: |
| 1          | right          | 1           | top-right   |
| 2          | top-right      | 2           | top         |
| 3          | top-left       | 3           | top-left    |
| 4          | left           | 4           | bottom-left |
| 5          | bottom-left    | 5           | bottom      |
| 6          | bottom-right   | 6           | bottom-right|

"""
function TriangularMesh(; dx=1e-9, dz=1e-9, nx=1, ny=1, pbc="open")
    nz = 1
    dz = 1e-9
    n_total = nx * ny * nz
    ngbs = zeros(Int32, 6, n_total)
    nngbs = zeros(Int32, 6, n_total)
    pbc_x = 'x' in pbc ? true : false
    pbc_y = 'y' in pbc ? true : false
    pbc_z = 'z' in pbc ? true : false
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i + 1, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #right
        ngbs[2, id] = indexpbc(i + 1, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #top_right
        ngbs[3, id] = indexpbc(i, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #top_left
        ngbs[4, id] = indexpbc(i - 1, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #left
        ngbs[5, id] = indexpbc(i - 1, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #bottom_left
        ngbs[6, id] = indexpbc(i, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #bottom_right

        nngbs[1, id] = indexpbc(i + 2, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #top_right
        nngbs[2, id] = indexpbc(i + 1, j + 2, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #top
        nngbs[3, id] = indexpbc(i - 1, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #top_left
        nngbs[4, id] = indexpbc(i - 2, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #bottom_left
        nngbs[5, id] = indexpbc(i - 1, j - 2, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #bottom
        nngbs[6, id] = indexpbc(i + 1, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)  #bottom_right
    end
    dy = dx * sqrt(3) / 2

    ngbs = kernel_array(ngbs)
    nngbs = kernel_array(nngbs)

    return TriangularMesh(dx, dy, dz, nx, ny, nz, n_total, 6, 6, ngbs, nngbs, pbc_x, pbc_y,
                          pbc_z)
end

mutable struct CubicMesh <: AtomisticMesh
    dx::Float64
    dy::Float64
    dz::Float64
    nx::Int64
    ny::Int64
    nz::Int64
    n_total::Int64
    n_ngbs::Int64  #number of the nearest neighbours 
    n_ngbs2::Int64  #number of next-nearest neighbours
    n_ngbs3::Int64  #number of next-next-nearest neighbours
    n_ngbs4::Int64  #number of next-next-next-nearest neighbours
    ngbs::AbstractArray{Int32,2}  # nearest neighbours
    ngbs2::AbstractArray{Int32,2}  # next-nearest neighbours 
    ngbs3::AbstractArray{Int32,2}  # next-next-nearest neighbours
    ngbs4::AbstractArray{Int32,2}  # next-next-next-nearest neighbours
    xperiodic::Bool
    yperiodic::Bool
    zperiodic::Bool
end

function Base.size(mesh::CubicMesh)
    return mesh.nx, mesh.ny, mesh.nz
end

function Base.length(mesh::Mesh)
    return mesh.n_total
end

@doc raw"""	 
	CubicMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")

Create a simple cubic mesh, in which each cell is only connected to six nearest neighbors 
and forming a simple cube structure.
"""
function CubicMesh(; dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
    n_total = nx * ny * nz
    ngbs = zeros(Int32, 6, n_total)

    pbc_x = 'x' in pbc ? true : false
    pbc_y = 'y' in pbc ? true : false
    pbc_z = 'z' in pbc ? true : false
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i - 1, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[2, id] = indexpbc(i + 1, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[3, id] = indexpbc(i, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[4, id] = indexpbc(i, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[5, id] = indexpbc(i, j, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[6, id] = indexpbc(i, j, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
    end

    ngbs = kernel_array(ngbs)
    empty_ngbs = Array{Int32}(undef, 1, 0)
    return CubicMesh(dx, dy, dz, nx, ny, nz, n_total, 6, 12, 8, 6, ngbs, empty_ngbs,
                     empty_ngbs, empty_ngbs, pbc_x, pbc_y, pbc_z)
end

"""
compute 2nd next-nearest neighbours for CubicMesh
"""
function compute_2nd_ngbs(mesh::CubicMesh)
    ngbs = zeros(Int32, 12, mesh.n_total)
    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz
    pbc_x = mesh.xperiodic
    pbc_y = mesh.yperiodic
    pbc_z = mesh.zperiodic
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i + 1, j, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[2, id] = indexpbc(i, j + 1, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[3, id] = indexpbc(i - 1, j, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[4, id] = indexpbc(i, j - 1, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[5, id] = indexpbc(i + 1, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[6, id] = indexpbc(i + 1, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[7, id] = indexpbc(i - 1, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[8, id] = indexpbc(i - 1, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[9, id] = indexpbc(i + 1, j, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[10, id] = indexpbc(i, j + 1, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[11, id] = indexpbc(i - 1, j, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[12, id] = indexpbc(i, j - 1, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
    end

    mesh.n_ngbs2 = 12
    mesh.ngbs2 = kernel_array(ngbs)
    return nothing
end

"""
compute 3rd next-next-nearest neighbours for CubicMesh
"""
function compute_3rd_ngbs(mesh::CubicMesh)
    ngbs = zeros(Int32, 8, mesh.n_total)
    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz
    pbc_x = mesh.xperiodic
    pbc_y = mesh.yperiodic
    pbc_z = mesh.zperiodic
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i + 1, j - 1, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[2, id] = indexpbc(i + 1, j + 1, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[3, id] = indexpbc(i - 1, j + 1, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[4, id] = indexpbc(i - 1, j - 1, k - 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[5, id] = indexpbc(i + 1, j - 1, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[6, id] = indexpbc(i + 1, j + 1, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[7, id] = indexpbc(i - 1, j + 1, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[8, id] = indexpbc(i - 1, j - 1, k + 1, nx, ny, nz, pbc_x, pbc_y, pbc_z)
    end

    mesh.n_ngbs3 = 8
    mesh.ngbs3 = kernel_array(ngbs)
    return nothing
end

"""
compute 4th next-next-next-nearest neighbours for CubicMesh
"""
function compute_4th_ngbs(mesh::CubicMesh)
    ngbs = zeros(Int32, 6, mesh.n_total)
    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz
    pbc_x = mesh.xperiodic
    pbc_y = mesh.yperiodic
    pbc_z = mesh.zperiodic
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i - 2, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[2, id] = indexpbc(i + 2, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[3, id] = indexpbc(i, j - 2, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[4, id] = indexpbc(i, j + 2, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[5, id] = indexpbc(i, j, k - 2, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[6, id] = indexpbc(i, j, k + 2, nx, ny, nz, pbc_x, pbc_y, pbc_z)
    end

    mesh.n_ngbs4 = 6
    mesh.ngbs4 = kernel_array(ngbs)
    return nothing
end

# special for 2D-lattice (nz == 1), where the next-next-nearest points are totally different from the 3D case. 
struct SquareMesh <: AtomisticMesh
    dx::Float64
    dy::Float64
    dz::Float64
    nx::Int64
    ny::Int64
    nz::Int64
    n_total::Int64
    n_ngbs::Int64
    n_ngbs2::Int64
    n_ngbs3::Int64
    ngbs::AbstractArray{Int32,2}  #  nearest neighbours
    ngbs2::AbstractArray{Int32,2}  # next-nearest neighbours(warnning: the corresponding interaction function is named as "next_exch", not "nnexch")
    ngbs3::AbstractArray{Int32,2}  # next-next-nearest neighbours
    xperiodic::Bool
    yperiodic::Bool
    zperiodic::Bool
end

"""
    SquareMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, pbc="open")

Create a SquareMesh which is two dimensional version of the CubicMesh. The 
next-next-nearest neighbours are totally different from the CubicMesh. 
"""
function SquareMesh(; dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
    ngbs = zeros(Int32, 4, nx * ny * nz)
    nngbs = zeros(Int32, 4, nx * ny * nz)
    nnngbs = zeros(Int32, 4, nx * ny * nz)
    pbc_x = 'x' in pbc ? true : false
    pbc_y = 'y' in pbc ? true : false
    pbc_z = 'z' in pbc ? true : false
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1, id] = indexpbc(i - 1, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[2, id] = indexpbc(i + 1, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[3, id] = indexpbc(i, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        ngbs[4, id] = indexpbc(i, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)

        nngbs[1, id] = indexpbc(i + 1, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        nngbs[2, id] = indexpbc(i + 1, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        nngbs[3, id] = indexpbc(i - 1, j + 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        nngbs[4, id] = indexpbc(i - 1, j - 1, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)

        nnngbs[1, id] = indexpbc(i + 2, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        nnngbs[2, id] = indexpbc(i, j + 2, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        nnngbs[3, id] = indexpbc(i - 2, j, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
        nnngbs[4, id] = indexpbc(i, j - 2, k, nx, ny, nz, pbc_x, pbc_y, pbc_z)
    end
    n_total = nx * ny * nz

    ngbs = kernel_array(ngbs)
    nngbs = kernel_array(nngbs)
    nnngbs = kernel_array(nnngbs)
    return CubicMeshGPU(dx, dy, dz, nx, ny, nz, n_total, 4, 4, 4, ngbs, nngbs, nnngbs,
                        pbc_x, pbc_y, pbc_z)
end

#We defined a cylindrical tube mesh along the +z direction
struct CylindricalTubeMesh <: AtomisticMesh
    dz::Float64
    R::Float64
    nz::Int64
    nr::Int64
    n_total::Int64
    n_ngbs::Int64
    ngbs::AbstractArray{Int32,2}  #nearest neighbours
    zperiodic::Bool
    coordinates::Array{Float64,2}  #coordinates array
end

"""
    CylindricalTubeMesh(;dz=1e-9, R=20e-9, nz=1, nr=10, pbc="open")

Create a cylindrical tube mesh along the +z direction. 

The spins are located on the cylindrical tube uniformly, and are indexed as follows:
```julia
  id = index(i, 1, k, nr, 1, nz)
```
which means that the spins are labelled in a ring firstly, then `nz` is the number of rings.  
The coordinates of the spins at each ring are given as `(R cos(2*pi*(i-1)/nr), R sin(2*pi*(nr-1)/nr))`.

The nearest neighbours are indexed as follows:

  |  1      2        3       4         
  |left   right    bottom   top |
"""
function CylindricalTubeMesh(; dz=1e-9, R=20e-9, nz=1, nr=10, pbc="open")
    ngbs = zeros(Int32, 4, nr * nz)
    coordinates = zeros(Float64, 3, nr * nz)
    pbc_z = 'z' in pbc ? true : false
    delta = 2 * pi / nr
    for k in 1:nz, i in 1:nr
        id = index(i, 1, k, nr, 1, nz)
        coordinates[1, id] = R * cos((i - 1) * delta)
        coordinates[2, id] = R * sin((i - 1) * delta)
        coordinates[3, id] = (k - 0.5) * dz
        ngbs[1, id] = indexpbc(i - 1, 1, k, nr, 1, nz, true, false, pbc_z)
        ngbs[2, id] = indexpbc(i + 1, 1, k, nr, 1, nz, true, false, pbc_z)
        ngbs[3, id] = indexpbc(i, 1, k - 1, nr, 1, nz, true, false, pbc_z)
        ngbs[4, id] = indexpbc(i, 1, k + 1, nr, 1, nz, true, false, pbc_z)
    end
    n_total = nr * nz
    n_ngbs = 4
    ngbs = kernel_array(ngbs)
    return CylindricalTubeMesh(dz, R, nz, nr, n_total, n_ngbs, ngbs, pbc_z, coordinates)
end

# We assume this mesh can be used for NiO
# We consider 4 Ni atoms in each cell, a cell has a size (dx, dy, dz)
# Four Ni atoms are located at (0,0,0), (dx/2, dy/2, 0), (0, dy/2, dz/2), (dx/2, 0, dz/2)
struct FccMesh <: AtomisticMesh
    dx::Float64
    dy::Float64
    dz::Float64
    nx::Int64
    ny::Int64
    nz::Int64
    n_total::Int64
    n_ngbs::Int64
    nn_ngbs::Int64
    ngbs::AbstractArray{Int32,2}  #  nearest neighbours
    nngbs::AbstractArray{Int32,2}  # next-nearest neighbours
    xperiodic::Bool
    yperiodic::Bool
    zperiodic::Bool
    coordinates::Array{Float64,2}  #coordinates array
end

function distance_x(x0, x1, L)
    x = x0 - x1
    if L <= 0
        return abs(x)
    end
    if x >= 0
        return min(x, abs(x - L))
    elseif x < 0
        return min(-x, L + x)
    end
end

function FccMesh(; dx=1e-9, dy=1e-9, dz=1e-9, nx=2, ny=2, nz=2, pbc="xyz")
    coordinates = zeros(Float64, 3, nx * ny * nz * 4)
    pbc_x = 'x' in pbc ? true : false
    pbc_y = 'y' in pbc ? true : false
    pbc_z = 'z' in pbc ? true : false

    for k in 1:nz, j in 1:ny, i in 1:nx
        id = 4 * (index(i, j, k, nx, ny, nz) - 1) + 1
        coordinates[1, id] = dx * i
        coordinates[2, id] = dy * j
        coordinates[3, id] = dz * k
        coordinates[1, id + 1] = dx * (i + 0.5)
        coordinates[2, id + 1] = dy * (j + 0.5)
        coordinates[3, id + 1] = dz * k
        coordinates[1, id + 2] = dx * i
        coordinates[2, id + 2] = dy * (j + 0.5)
        coordinates[3, id + 2] = dz * (k + 0.5)
        coordinates[1, id + 3] = dx * (i + 0.5)
        coordinates[2, id + 3] = dy * j
        coordinates[3, id + 3] = dz * (k + 0.5)
    end

    N = nx * ny * nz * 4
    cds = coordinates
    Lx = pbc_x ? nx * dx : -1
    Ly = pbc_y ? ny * dy : -1
    Lz = pbc_z ? nz * dz : -1
    d0 = sqrt(dx^2 + dy^2 + dz^2) / 2
    d1 = (dx + dy + dz) / 3.0
    ngbs_array = []
    nngbs_array = []
    for i in 1:N
        _ngbs = Int32[]
        _nngbs = Int32[]
        for j in 1:N
            if j == i
                continue
            end

            x = distance_x(cds[1, i], cds[1, j], Lx)
            y = distance_x(cds[2, i], cds[2, j], Ly)
            z = distance_x(cds[3, i], cds[3, j], Lz)

            d = sqrt(x^2 + y^2 + z^2)
            if d < d0
                push!(_ngbs, j)
            elseif d < d1 * 1.1
                push!(_nngbs, j)
            end
        end
        push!(ngbs_array, _ngbs)
        push!(nngbs_array, _nngbs)
    end

    n_ngbs = length(ngbs_array[1])
    n_nngbs = length(nngbs_array[1])
    ngbs = zeros(Int32, n_ngbs, N)
    nngbs = zeros(Int32, n_nngbs, N)
    for i in 1:N
        ngbs[:, i] .= ngbs_array[i]
        nngbs[:, i] .= nngbs_array[i]
    end
    ngbs = kernel_array(ngbs)
    nngbs = kernel_array(nngbs)
    return FccMesh(dx, dy, dz, nx, ny, nz, N, n_ngbs, n_nngbs, ngbs, nngbs, pbc_x, pbc_y,
                   pbc_z, cds)
end
