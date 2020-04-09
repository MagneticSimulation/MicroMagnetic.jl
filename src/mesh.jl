abstract type Mesh end

struct FDMesh <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  volume::Float64
  ngbs::Array{Int64, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end


struct CubicMesh <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  volume::Float64
  ngbs::Array{Int64, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

struct TriangularMesh <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  ngbs::Array{Int64, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

"""
    FDMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")

Create a FDMesh for given parameters. `pbc` could be any combination of "x", "y" and "z".
"""
function FDMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
    ngbs = zeros(Int64, 6, nx*ny*nz)
    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false
    zperiodic = 'z' in pbc ? true : false

    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1,id] = indexpbc(i-1,j,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[2,id] = indexpbc(i+1,j,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[3,id] = indexpbc(i,j-1,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[4,id] = indexpbc(i,j+1,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[5,id] = indexpbc(i,j,k-1, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[6,id] = indexpbc(i,j,k+1, nx,ny,nz, xperiodic, yperiodic, zperiodic)
    end
    volume = dx*dy*dz
    nxyz = nx*ny*nz
    return FDMesh(dx, dy, dz, nx, ny, nz, nxyz, volume, ngbs, xperiodic, yperiodic, zperiodic)
end

"""
    CubicMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
"""
function CubicMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
    ngbs = zeros(Int64, 6, nx*ny*nz)
    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false
    zperiodic = 'z' in pbc ? true : false

    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        ngbs[1,id] = indexpbc(i-1,j,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[2,id] = indexpbc(i+1,j,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[3,id] = indexpbc(i,j-1,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[4,id] = indexpbc(i,j+1,k, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[5,id] = indexpbc(i,j,k-1, nx,ny,nz, xperiodic, yperiodic, zperiodic)
        ngbs[6,id] = indexpbc(i,j,k+1, nx,ny,nz, xperiodic, yperiodic, zperiodic)
    end
    volume = dx*dy*dz
    nxyz = nx*ny*nz
    return CubicMesh(dx, dy, dz, nx, ny, nz, nxyz, volume, ngbs, xperiodic, yperiodic, zperiodic)
end

"""
Create a 3d triangular mesh.

The neighbours are indexed as counterclockwise of the given spin:

  |  1      2         3       4         5          6           7      8   |
  |right top_right top_left  left  bottom_left bottom_right   below  above |

"""
function TriangularMesh(;dx=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
    ngbs = zeros(Int64, 8, nx*ny*nz)
    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false
    zperiodic = 'z' in pbc ? true : false

    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i,j, k, nx, ny, nz)
        ngbs[1,id] = indexpbc(i+1, j, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #right
        ngbs[2,id] = indexpbc(i+1, j+1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #top_right
        ngbs[3,id] = indexpbc(i, j+1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #top_left
        ngbs[4,id] = indexpbc(i-1, j, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #left
        ngbs[5,id] = indexpbc(i-1, j-1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #bottom_left
        ngbs[6,id] = indexpbc(i, j-1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #bottom_right
        ngbs[7,id] = indexpbc(i, j, k-1, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #below
        ngbs[8,id] = indexpbc(i, j, k+1, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #above
    end
    dy  = dx*sqrt(3)/2
    return TriangularMesh(dx, dy, dz, nx, ny, nz, nx*ny*nz, ngbs, xperiodic, yperiodic, zperiodic)
end


@inline function index(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64)
  if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
    return -1
  end
  return (k-1) * nx*ny + (j-1) * nx + i
end

@inline function _x_plus_one(i::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, xperiodic::Bool)::Int64
    if i < nx || xperiodic
        return (i==nx) ? index +1 - nx : index +1
    end
    return -1
end

@inline function _x_minus_one(i::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, xperiodic::Bool)::Int64
    if i > 1 || xperiodic
        return (i==1) ? index - 1 + nx : index -1
    end
    return -1
end

@inline function _y_plus_one(j::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, yperiodic::Bool)::Int64
    if j < ny || yperiodic
        return (j==ny) ? index + nx - nx*ny : index + nx
    end
    return -1
end

@inline function _y_minus_one(j::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, yperiodic::Bool)::Int64
    if j > 1 || yperiodic
        return (j==1) ? index - nx + nx*ny : index - nx
    end
    return -1
end

@inline function _z_plus_one(k::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, zperiodic::Bool)::Int64
    if k<nz || zperiodic
        return (k==nz) ? index + nx*ny*(1-nz) : index + nx*ny
	end
    return -1
end

@inline function _z_minus_one(k::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, zperiodic::Bool)::Int64
    if k>1 || zperiodic
        return (k==1) ? index + nx*ny*(nz-1) : index - nx*ny
	end
    return -1
end


function indexpbc(i::Int64, j::Int64, k::Int64,
                  nx::Int64, ny::Int64, nz::Int64,
                  xperiodic::Bool, yperiodic::Bool, zperiodic::Bool)
  if xperiodic
    if i < 1
      i += nx
    elseif i > nx
      i -= nx
    end
  end

  if yperiodic
    if j < 1
      j += ny
    elseif j > ny
      j -= ny
    end
  end

  if zperiodic
    if k < 1
      k += nz
    elseif k > nz
      k -= nz
    end
  end
  return index(i,j,k, nx, ny, nz)
end
