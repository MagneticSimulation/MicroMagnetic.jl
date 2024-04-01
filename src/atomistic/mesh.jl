abstract type AtomicMeshGPU <: MeshGPU end

struct TriangularMeshGPU <: AtomicMeshGPU
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  n_total::Int64
  n_ngbs::Int64
  ngbs::CuArray{Int32, 2}
  nngbs::CuArray{Int32, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

@doc raw"""
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


struct CubicMeshGPU <: AtomicMeshGPU
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  n_total::Int64
  n_ngbs::Int64
  nn_ngbs::Int64
  nnn_ngbs::Int64
  ngbs::CuArray{Int32, 2}  #  nearest neighbours
  nngbs::CuArray{Int32, 2}  # next-nearest neighbours(warnning: the corresponding interaction function is named as "next_exch", not "nnexch")
  nnngbs::CuArray{Int32, 2}  # next-next-nearest neighbours
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

@doc raw"""
Create a CubicMesh 
	 
	CubicMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")

"""
function CubicMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int32,6,nx*ny*nz)
  nngbs = zeros(Int32,12,nx*ny*nz)
  nnngbs = zeros(Int32, 8, nx*ny*nz)
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

    nngbs[1,id] = indexpbc(i+1,j,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[2,id] = indexpbc(i,j+1,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[3,id] = indexpbc(i-1,j,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[4,id] = indexpbc(i,j-1,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[5,id] = indexpbc(i+1,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[6,id] = indexpbc(i+1,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[7,id] = indexpbc(i-1,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[8,id] = indexpbc(i-1,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[9,id] = indexpbc(i+1,j,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[10,id] = indexpbc(i,j+1,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[11,id] = indexpbc(i-1,j,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[12,id] = indexpbc(i,j-1,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)

    nnngbs[1,id] = indexpbc(i+1,j-1,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[2,id] = indexpbc(i+1,j+1,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[3,id] = indexpbc(i-1,j+1,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[4,id] = indexpbc(i-1,j-1,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[5,id] = indexpbc(i+1,j-1,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[6,id] = indexpbc(i+1,j+1,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[7,id] = indexpbc(i-1,j+1,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[8,id] = indexpbc(i-1,j-1,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
  end
  n_total = nx*ny*nz
  return CubicMeshGPU(dx, dy, dz, nx, ny, nz, n_total, 6, 12, 8, CuArray(ngbs), CuArray(nngbs), CuArray(nnngbs), xperiodic, yperiodic, zperiodic)
end

struct CubicMesh4NGPU <: AtomicMeshGPU
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  n_total::Int64
  n_ngbs::Int64
  nn_ngbs::Int64
  nnn_ngbs::Int64
  nnnn_ngbs::Int64
  ngbs::CuArray{Int32, 2}  #  1st nearest neighbours
  nngbs::CuArray{Int32, 2}  # 2nd next-nearest neighbours(warnning: the corresponding interaction function is named as "next_exch", not "nnexch")
  nnngbs::CuArray{Int32,2}  # 3rd next-next-nearest neighbours
  nnnngbs::CuArray{Int32,2}  # 4th next-next-next-nearest neighbours
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

function CubicMesh4NGPU(; dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int32, 6, nx * ny * nz)
  nngbs = zeros(Int32, 12, nx * ny * nz)
  nnngbs = zeros(Int32, 8, nx * ny * nz)
  nnnngbs = zeros(Int32, 6, nx * ny * nz)
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

    nngbs[1, id] = indexpbc(i + 1, j, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[2, id] = indexpbc(i, j + 1, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[3, id] = indexpbc(i - 1, j, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[4, id] = indexpbc(i, j - 1, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[5, id] = indexpbc(i + 1, j - 1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[6, id] = indexpbc(i + 1, j + 1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[7, id] = indexpbc(i - 1, j + 1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[8, id] = indexpbc(i - 1, j - 1, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[9, id] = indexpbc(i + 1, j, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[10, id] = indexpbc(i, j + 1, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[11, id] = indexpbc(i - 1, j, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nngbs[12, id] = indexpbc(i, j - 1, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)

    nnngbs[1, id] = indexpbc(i + 1, j - 1, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[2, id] = indexpbc(i + 1, j + 1, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[3, id] = indexpbc(i - 1, j + 1, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[4, id] = indexpbc(i - 1, j - 1, k - 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[5, id] = indexpbc(i + 1, j - 1, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[6, id] = indexpbc(i + 1, j + 1, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[7, id] = indexpbc(i - 1, j + 1, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnngbs[8, id] = indexpbc(i - 1, j - 1, k + 1, nx, ny, nz, xperiodic, yperiodic, zperiodic)

    nnnngbs[1, id] = indexpbc(i - 2, j, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnnngbs[2, id] = indexpbc(i + 2, j, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnnngbs[3, id] = indexpbc(i, j - 2, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnnngbs[4, id] = indexpbc(i, j + 2, k, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnnngbs[5, id] = indexpbc(i, j, k - 2, nx, ny, nz, xperiodic, yperiodic, zperiodic)
    nnnngbs[6, id] = indexpbc(i, j, k + 2, nx, ny, nz, xperiodic, yperiodic, zperiodic)
  end
  n_total = nx * ny * nz
  return CubicMesh4NGPU(dx, dy, dz, nx, ny, nz, n_total, 6, 12, 8, 6, CuArray(ngbs), CuArray(nngbs), CuArray(nnngbs), CuArray(nnnngbs), xperiodic, yperiodic, zperiodic)
end


#   special for 2D-lattice(nz == 1), where the next-next-nearest points are totally different from the 3D case. 
struct SquareMeshGPU <: AtomicMeshGPU
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  n_total::Int64
  n_ngbs::Int64
  nn_ngbs::Int64
  nnn_ngbs::Int64
  ngbs::CuArray{Int32, 2}  #  nearest neighbours
  nngbs::CuArray{Int32, 2}  # next-nearest neighbours(warnning: the corresponding interaction function is named as "next_exch", not "nnexch")
  nnngbs::CuArray{Int32, 2}  #    next-next-nearest neighbours
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end


function SquareMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int32,4,nx*ny*nz)
  nngbs = zeros(Int32,4,nx*ny*nz)
  nnngbs = zeros(Int32, 4, nx*ny*nz)
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false
  for k = 1:nz, j = 1:ny, i=1:nx
    id = index(i,j,k, nx, ny, nz)
    ngbs[1,id] = indexpbc(i-1,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[2,id] = indexpbc(i+1,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[3,id] = indexpbc(i,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[4,id] = indexpbc(i,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)

    nngbs[1,id] = indexpbc(i+1,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[2,id] = indexpbc(i+1,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[3,id] = indexpbc(i-1,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nngbs[4,id] = indexpbc(i-1,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)

    nnngbs[1,id] = indexpbc(i+2,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[2,id] = indexpbc(i,j+2,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[3,id] = indexpbc(i-2,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    nnngbs[4,id] = indexpbc(i,j-2,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
  end
  n_total = nx*ny*nz
  return CubicMeshGPU(dx, dy, dz, nx, ny, nz, n_total, 4, 4, 4, CuArray(ngbs), CuArray(nngbs), CuArray(nnngbs), xperiodic, yperiodic, zperiodic)
end

#We defined a cylindrical tube mesh along the +z direction
struct CylindricalTubeMeshGPU <: AtomicMeshGPU
  dz::Float64
  R::Float64
  nz::Int64
  nr::Int64
  n_total::Int64
  n_ngbs::Int64
  ngbs::CuArray{Int32, 2}  #nearest neighbours
  zperiodic::Bool
  coordinates::Array{Float64, 2}  #coordinates array
end

"""
    CylindricalTubeMeshGPU(;dz=1e-9, R=20e-9, nz=1, nr=10, pbc="open")

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
function CylindricalTubeMeshGPU(;dz=1e-9, R=20e-9, nz=1, nr=10, pbc="open")
  ngbs = zeros(Int32, 4, nr*nz)
  coordinates = zeros(Float64, 3, nr*nz)
  zperiodic = 'z' in pbc ? true : false
  delta = 2*pi/nr
  for k = 1:nz, i=1:nr
    id = index(i, 1, k, nr, 1, nz)
    coordinates[1, id] = R*cos((i-1)*delta) 
    coordinates[2, id] = R*sin((i-1)*delta) 
    coordinates[3, id] = (k-0.5)*dz
    ngbs[1,id] = indexpbc(i-1,1,k,nr,1,nz, true, false, zperiodic)
    ngbs[2,id] = indexpbc(i+1,1,k,nr,1,nz, true, false, zperiodic)
    ngbs[3,id] = indexpbc(i,1,k-1,nr,1,nz, true, false, zperiodic)
    ngbs[4,id] = indexpbc(i,1,k+1,nr,1,nz, true, false, zperiodic)
  end
  n_total = nr*nz
  n_ngbs = 4
  return CylindricalTubeMeshGPU(dz, R, nz, nr, n_total, n_ngbs, CuArray(ngbs), zperiodic, coordinates)
end


# We assume this mesh can be used for NiO
# We consider 4 Ni atoms in each cell, a cell has a size (dx, dy, dz)
# Four Ni atoms are located at (0,0,0), (dx/2, dy/2, 0), (0, dy/2, dz/2), (dx/2, 0, dz/2)
struct FccMeshGPU <: AtomicMeshGPU
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  n_total::Int64  # it should be n_spins, we will change it later.
  n_ngbs::Int64
  nn_ngbs::Int64
  ngbs::CuArray{Int32, 2}  #  nearest neighbours
  nngbs::CuArray{Int32, 2} 
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
  coordinates::Array{Float64, 2}  #coordinates array
end

function distance_x(x0, x1, L)
  x = x0 - x1
  if L <= 0
    return abs(x)
  end
  if x>=0
    return min(x, abs(x-L))
  elseif x<0
    return min(-x, L+x)
  end
end

function FccMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=2, ny=2, nz=2, pbc="xyz")
  
  coordinates = zeros(Float64, 3, nx*ny*nz*4)
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false
  
  for k = 1:nz, j=1:ny, i=1:nx
    id = 4*(index(i, j, k, nx, ny, nz)-1)+1
    coordinates[1, id] = dx*i
    coordinates[2, id] = dy*j
    coordinates[3, id] = dz*k
    coordinates[1, id+1] = dx*(i+0.5)
    coordinates[2, id+1] = dy*(j+0.5)
    coordinates[3, id+1] = dz*k
    coordinates[1, id+2] = dx*i
    coordinates[2, id+2] = dy*(j+0.5)
    coordinates[3, id+2] = dz*(k+0.5)
    coordinates[1, id+3] = dx*(i+0.5)
    coordinates[2, id+3] = dy*j
    coordinates[3, id+3] = dz*(k+0.5)
  end

  N = nx*ny*nz*4
  cds = coordinates
  Lx = xperiodic ? nx*dx : -1
  Ly = yperiodic ? ny*dy : -1
  Lz = zperiodic ? nz*dz : -1
  d0 = sqrt(dx^2+dy^2+dz^2)/2
  d1 = (dx+dy+dz)/3.0
  ngbs_array = []
  nngbs_array = []
  for i = 1:N
    _ngbs = Int32[]
    _nngbs = Int32[]
    for j=1:N
      if j == i
        continue 
      end

      x = distance_x(cds[1, i], cds[1, j], Lx)
      y = distance_x(cds[2, i], cds[2, j], Ly)
      z = distance_x(cds[3, i], cds[3, j], Lz)
      
      d = sqrt(x^2+y^2+z^2)
      if d<d0
        push!(_ngbs, j)
      elseif d < d1*1.1
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
  for i = 1:N
    ngbs[:, i] .= ngbs_array[i]
    nngbs[:, i] .= nngbs_array[i]
  end
  return FccMeshGPU(dx, dy, dz, nx, ny, nz, N, n_ngbs, n_nngbs, cu(ngbs), cu(nngbs), xperiodic, yperiodic, zperiodic, cds)
end