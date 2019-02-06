struct Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  ngbs::Array{Int64}
end

function index(i, j, k, nx, ny, nz)
  if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
    return -1
  end
  return (k-1) * nx*ny + (j-1) * nx + i
end

function create_mesh(;dx=1.0, dy=1.0, dz=1.0, nx=1, ny=1, nz=1)
  ngbs = zeros(Int64,6,nx*ny*nz)
  for i = 1:nx
    for j = 1:ny
      for k = 1:nz
        id = (k-1) * nx*ny + (j-1) * nx + i
        ngbs[1,id] = index(i-1,j,k,nx,ny,nz)
        ngbs[2,id] = index(i+1,j,k,nx,ny,nz)
        ngbs[3,id] = index(i,j-1,k,nx,ny,nz)
        ngbs[4,id] = index(i,j+1,k,nx,ny,nz)
        ngbs[5,id] = index(i,j,k-1,nx,ny,nz)
        ngbs[6,id] = index(i,j,k+1,nx,ny,nz)
      end
    end
  end
  return Mesh(dx, dy, dz, nx, ny, nz, ngbs)
end
