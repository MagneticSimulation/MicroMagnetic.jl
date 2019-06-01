struct TriangularMesh <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  ngbs::CuArray{Int32, 2}
  xperiodic::Bool
  yperiodic::Bool
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
Create a 2d triangular mesh.

The neighbours are indexed as counterclockwise of the given spin:

  |  1      2         3       4         5          6        |
  |right top_right top_left  left  bottom_left bottom_right |

"""
function TriangularMesh(;dx=1e-9, nx=3, ny=2, pbc="open")

    ngbs = zeros(Int32, 6, nx*ny)
    xperiodic = 'x' in pbc ? true : false
    yperiodic = 'y' in pbc ? true : false

    for j = 1:ny, i=1:nx
        id = index(i, j, nx, ny)
        ngbs[1,id] = indexpbc(i+1, j, nx, ny, xperiodic, yperiodic)  #right
        ngbs[2,id] = indexpbc(i, j+1, nx, ny, xperiodic, yperiodic)  #top_right
        ngbs[3,id] = indexpbc(i-1, j+1, nx, ny, xperiodic, yperiodic)  #top_left
        ngbs[4,id] = indexpbc(i-1, j, nx, ny, xperiodic, yperiodic)  #left
        ngbs[5,id] = indexpbc(i, j-1, nx, ny, xperiodic, yperiodic)  #bottom_left
        ngbs[6,id] = indexpbc(i+1, j-1, nx, ny, xperiodic, yperiodic)  #bottom_right
    end
    dy  = dx*sqrt(3)/2
    return TriangularMesh(dx, dy, dx, nx, ny, 1, nx*ny, CuArray(ngbs), xperiodic, yperiodic)
end
