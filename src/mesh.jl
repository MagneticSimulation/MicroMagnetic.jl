@inline function index(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64)
  if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
    return -1
  end
  return (k-1) * nx*ny + (j-1) * nx + i
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

function FDMesh(;dx=1.0, dy=1.0, dz=1.0, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int64,6,nx*ny*nz)
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false

  for k = 1:nz, j = 1:ny, i=1:nx
    id = index(i, j, k, nx, ny, nz)
    ngbs[1,id] = indexpbc(i-1,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[2,id] = indexpbc(i+1,j,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[3,id] = indexpbc(i,j-1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[4,id] = indexpbc(i,j+1,k,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[5,id] = indexpbc(i,j,k-1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
    ngbs[6,id] = indexpbc(i,j,k+1,nx,ny,nz, xperiodic, yperiodic, zperiodic)
  end
  volume = dx*dy*dz
  nxyz = nx*ny*nz
  return FDMesh(dx, dy, dz, nx, ny, nz, nxyz, volume, ngbs, xperiodic, yperiodic, zperiodic)
end

function CubicMesh(;a=1.0, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int64,6,nx*ny*nz)
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
  end
  nxyz = nx*ny*nz
  return CubicMesh(a, nx, ny, nz, nxyz, ngbs, xperiodic, yperiodic, zperiodic)
end
