abstract type Mesh end
abstract type MeshGPU <: Mesh end
abstract type AtomicMeshGPU <: MeshGPU end

struct TriMesh3DGPU <: AtomicMeshGPU
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  TotRdN::Int64 #有多少近邻 2:xy,z
  n_ngbs::CuArray{Int32, 1} #第几近邻到第几个，0,6,8
  ngbs::CuArray{Int32,2}  #按照123近邻顺序写入
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

"""
Create a 3d triangular mesh.
a and b angel 60
The nearest neighbours are indexed as counterclockwise of the given spin:

  |  1      2         3       4         5          6           7      8   |
  |right top_right top_left  left  bottom_left bottom_right   below  above |



"""
function TriMesh3DGPU(;dx=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="xyz")
    TotRdN=2  #xy,z
    nn = zeros(Int32, 3)
    nn[1]=Int32(0)
    nn[2]=Int32(6)
    nn[3]=Int32(8)
    # nn[4]=Int32(18)
    ngbs = zeros(Int32, 8, nx*ny*nz)
    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false
    zperiodic = 'z' in pbc ? true : false
    for c = 1:nz, b = 1:ny, a=1:nx
        id = index(a,b, c, nx, ny, nz)
        ngbs[1,id] = indexpbc(a+1,b,  c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #right
        ngbs[2,id] = indexpbc(a+1,b+1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top_right
        ngbs[3,id] = indexpbc(a  ,b+1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top_left
        ngbs[4,id] = indexpbc(a-1,b,  c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #left
        ngbs[5,id] = indexpbc(a-1,b-1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom_left
        ngbs[6,id] = indexpbc(a  ,b-1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom_right

        ngbs[7,id] = indexpbc(a  ,b  ,c-1, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #below
        ngbs[8,id] = indexpbc(a  ,b  ,c+1, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #above

        # ngbs[7, id] = indexpbc(a+2,b+1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top_right
        # ngbs[8, id] = indexpbc(a+1,b+2,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top
        # ngbs[9, id] = indexpbc(a-1,b+1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top_left
        # ngbs[10,id] = indexpbc(a-2,b-1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom_left
        # ngbs[11,id] = indexpbc(a-1,b-2,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom
        # ngbs[12,id] = indexpbc(a+1,b-1,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom_right
        
        # ngbs[13,id] = indexpbc(a+2,b,  c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #right
        # ngbs[14,id] = indexpbc(a+2,b+2,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top_right
        # ngbs[15,id] = indexpbc(a  ,b+2,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #top_left
        # ngbs[16,id] = indexpbc(a-2,b,  c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #left
        # ngbs[17,id] = indexpbc(a-2,b-2,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom_left
        # ngbs[18,id] = indexpbc(a  ,b-2,c, nx,ny,nz, xperiodic, yperiodic, zperiodic)  #bottom_right

    end
    dy  = dx*sqrt(3)/2
    return TriMesh3DGPU(dx, dy, dz, nx, ny, nz, nx*ny*nz, TotRdN,CuArray(nn), CuArray(ngbs), xperiodic, yperiodic, zperiodic)
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

@inline function index(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64)
  if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
    return -1
  end
  return (k-1) * nx*ny + (j-1) * nx + i
end