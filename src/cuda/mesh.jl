function FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc="open")
  nxyz = nx*ny*nz
  volume = dx*dy*dz
  Float = _cuda_using_double.x ? Float64 : Float32
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false
  return FDMeshGPU(Float(dx), Float(dy), Float(dz), nx, ny, nz, nxyz, xperiodic, yperiodic, zperiodic, Float(volume))
end

struct TriangularMeshGPU <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  n_ngbs::Int64
  ngbs::CuArray{Int32, 2}
  nngbs::CuArray{Int32, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

@inline function index(i::Int64, j::Int64, nx::Int64, ny::Int64)::Int32
  if i < 1 || j < 1 || j > ny || i > nx
    return -1
  end
  return (j-1) * nx + i
end

function indexpbc(i::Int64, j::Int64, nx::Int64, ny::Int64, xperiodic::Bool, yperiodic::Bool)::Int32
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
  return index(i, j, nx, ny)
end

"""
Create a 3d triangular mesh.

The nearest neighbours are indexed as counterclockwise of the given spin:

  |  1      2         3       4         5          6           7      8   |
  |right top_right top_left  left  bottom_left bottom_right   below  above |

and the next-nearest neighbours are

  |  1       2        3       4         5          6           7      8   |
  |top_right top  top_left  bottom_left bottom bottom_right   below  above |

"""
function TriangularMeshGPU(;dx=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
    ngbs = zeros(Int32, 8, nx*ny*nz)
    nngbs = zeros(Int32, 8, nx*ny*nz)
    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false
    zperiodic = 'z' in pbc ? true : false
#modi by Kong
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

        nngbs[1,id] = indexpbc(i+2, j+1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #top_right
        nngbs[2,id] = indexpbc(i+1, j+2, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #top
        nngbs[3,id] = indexpbc(i-1, j+1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #top_left
        nngbs[4,id] = indexpbc(i-2, j-1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #bottom_left
        nngbs[5,id] = indexpbc(i-1, j-2, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #bottom
        nngbs[6,id] = indexpbc(i+1, j-1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #bottom_right
        nngbs[7,id] = indexpbc(i, j, k-2, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #below
        nngbs[8,id] = indexpbc(i, j, k+2, nx, ny, nz, xperiodic, yperiodic, zperiodic)  #above
    end
    dy  = dx*sqrt(3)/2
    return TriangularMeshGPU(dx, dy, dz, nx, ny, nz, nx*ny*nz, 8, CuArray(ngbs), CuArray(nngbs), xperiodic, yperiodic, zperiodic)
end


struct CubicMeshGPU <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  n_ngbs::Int64
  ngbs::CuArray{Int32, 2}
  nngbs::CuArray{Int32, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

function CubicMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int32,6,nx*ny*nz)
  nngbs = zeros(Int32,6,nx*ny*nz)
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false
  for k = 1:nz, j = 1:ny, i=1:nx
    id = index(i,j,k, nx, ny, nz)
    ngbs[1,id] = indexpbc(i-1,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[2,id] = indexpbc(i+1,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[3,id] = indexpbc(i,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[4,id] = indexpbc(i,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[5,id] = indexpbc(i,j,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[6,id] = indexpbc(i,j,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)

    nngbs[1,id] = indexpbc(i-2,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[2,id] = indexpbc(i+2,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[3,id] = indexpbc(i,j-2,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[4,id] = indexpbc(i,j+2,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[5,id] = indexpbc(i,j,k-2,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[6,id] = indexpbc(i,j,k+2,nx,ny,nz, xperiodic, yperiodic, zperiodic)
  end
  nxyz = nx*ny*nz
  return CubicMeshGPU(dx, dy, dz, nx, ny, nz, nxyz, 6, CuArray(ngbs), CuArray(nngbs), xperiodic, yperiodic, zperiodic)
end
