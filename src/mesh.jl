function index(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64)
  if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
    return -1
  end
  return (k-1) * nx*ny + (j-1) * nx + i
end

function indexpbc(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64, pbc::String)
  if 'x' in pbc
    if i < 1
      i += nx
    elseif i > nx
      i -= nx
    end
  end

  if 'y' in pbc
    if j < 1
      j += ny
    elseif j > ny
      j -= ny
    end
  end

  if 'z' in pbc
    if k < 1
      k += nz
    elseif k > nz
      k -= nz
    end
  end
  return index(i,j,k, nx, ny, nz)
end

function FDMesh(;dx=1.0, dy=1.0, dz=1.0, nx=1, ny=1, nz=1, unit_length=1.0, pbc="open")
  ngbs = zeros(Int64,6,nx*ny*nz)
  for k = 1:nz, j = 1:ny, i=1:nx
    id = index(i,j,k, nx, ny, nz)
    ngbs[1,id] = indexpbc(i-1,j,k,nx,ny,nz,pbc)
    ngbs[2,id] = indexpbc(i+1,j,k,nx,ny,nz,pbc)
    ngbs[3,id] = indexpbc(i,j-1,k,nx,ny,nz,pbc)
    ngbs[4,id] = indexpbc(i,j+1,k,nx,ny,nz,pbc)
    ngbs[5,id] = indexpbc(i,j,k-1,nx,ny,nz,pbc)
    ngbs[6,id] = indexpbc(i,j,k+1,nx,ny,nz,pbc)
  end
  volume = dx*dy*dz*unit_length^3
  nxyz = nx*ny*nz
  return FDMesh(dx, dy, dz, nx, ny, nz, nxyz, unit_length, volume, ngbs, pbc)
end

function CubicMesh(;a=1.0, nx=1, ny=1, nz=1, pbc="open")
  ngbs = zeros(Int64,6,nx*ny*nz)
  for k = 1:nz, j = 1:ny, i=1:nx
    id = index(i,j,k, nx, ny, nz)
    ngbs[1,id] = indexpbc(i-1,j,k,nx,ny,nz,pbc)
    ngbs[2,id] = indexpbc(i+1,j,k,nx,ny,nz,pbc)
    ngbs[3,id] = indexpbc(i,j-1,k,nx,ny,nz,pbc)
    ngbs[4,id] = indexpbc(i,j+1,k,nx,ny,nz,pbc)
    ngbs[5,id] = indexpbc(i,j,k-1,nx,ny,nz,pbc)
    ngbs[6,id] = indexpbc(i,j,k+1,nx,ny,nz,pbc)
  end
  nxyz = nx*ny*nz
  return CubicMesh(a, nx, ny, nz, nxyz, ngbs, pbc)
end
