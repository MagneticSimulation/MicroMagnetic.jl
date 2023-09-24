@inline function _x_plus_one(i::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, xperiodic::Bool,
                             Ms::CuDeviceArray{T, 1})::Int64 where {T<:AbstractFloat}
    if i < nx || xperiodic
        id = (i==nx) ? index +1 - nx : index +1
        if Ms[id]>0
          return id
        end
    end
    return -1
end

@inline function _x_minus_one(i::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, xperiodic::Bool,
                             Ms::CuDeviceArray{T, 1})::Int64 where {T<:AbstractFloat}
    if i > 1 || xperiodic
        id =  (i==1) ? index - 1 + nx : index -1
        if Ms[id]>0
          return id
        end
    end
    return -1
end

@inline function _y_plus_one(j::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, yperiodic::Bool,
                             Ms::CuDeviceArray{T, 1})::Int64 where {T<:AbstractFloat}
    if j < ny || yperiodic
        id = (j==ny) ? index + nx - nx*ny : index + nx
        if Ms[id]>0
          return id
        end
    end
    return -1
end

@inline function _y_minus_one(j::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, yperiodic::Bool,
                             Ms::CuDeviceArray{T, 1})::Int64 where {T<:AbstractFloat}
    if j > 1 || yperiodic
        id = (j==1) ? index - nx + nx*ny : index - nx
        if Ms[id]>0
          return id
        end
    end
    return -1
end

@inline function _z_plus_one(k::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, zperiodic::Bool,
                             Ms::CuDeviceArray{T, 1})::Int64 where {T<:AbstractFloat}
    if k<nz || zperiodic
        id = (k==nz) ? index + nx*ny*(1-nz) : index + nx*ny
        if Ms[id]>0
          return id
        end
	end
    return -1
end

@inline function _z_minus_one(k::Int64, index::Int64, nx::Int64, ny::Int64, nz::Int64, zperiodic::Bool,
                             Ms::CuDeviceArray{T, 1})::Int64 where {T<:AbstractFloat}
    if k>1 || zperiodic
        id = (k==1) ? index + nx*ny*(nz-1) : index - nx*ny
        if Ms[id]>0
          return id
        end
	end
    return -1
end

#compute a = a + b, which is slightly faster than a .+= b
function addto(a::CuArray{T,1}, b::CuArray{T,1})  where {T<:AbstractFloat}
    function __kernel!(a, b, n)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
      	if 0 < i <= n
           @inbounds a[i] += b[i]
        end
        return nothing
    end
    blk, thr = cudims(a)
    @cuda blocks=blk threads=thr __kernel!(a,b, length(a))
    return nothing
end

# compute the abs of a cuarray
function abs!(a::CuArray{T,1}, b::CuArray{T,1})  where {T<:AbstractFloat}
    function __kernel!(a, b, n)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        if 0<i<=n
            @inbounds a[i] = CUDA.abs(b[i])
        end
        return nothing
    end
    N = length(a)
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr __kernel!(a, b, N)
    return nothing
end

function abs!(a::CuArray{T,1})  where {T<:AbstractFloat}
    function __kernel!(a, n)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        if 0 < i <= n
            @inbounds a[i] = CUDA.abs(a[i])
        end
        return nothing
    end
    N = length(a)
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr __kernel!(a, N)
    return nothing
end

function  init_scalar!(v::CuArray{T, 1}, mesh::Mesh, init::Number) where {T<:AbstractFloat}
    v[:] .= init
    return true
end

function init_scalar!(v::CuArray{T, 1}, mesh::Mesh, init_fun::Function) where {T<:AbstractFloat}
    init_v = zeros(T, mesh.nxyz)
    init_scalar!(init_v, mesh, init_fun)
    copyto!(v, init_v)
    return true
end

function  init_scalar!(v::CuArray{T1, 1}, mesh::Mesh, init::Array{T2, 1}) where {T1,T2<:AbstractFloat}
    init_v  = zeros(T1, mesh.nxyz)
    init_v[:] .= init
    copyto!(v, init_v)
    return true
end

function init_vector!(v::Array{T, 1}, mesh::TriangularMeshGPU, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  dx,dy = mesh.dx, mesh.dy
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx, j = 1:mesh.ny, k = 1:mesh.nz
    id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
    vec_value = init(i,j,k, mesh.nx, mesh.ny, mesh.nz)
    if vec_value != nothing
      b[:, id] .= vec_value[:]
    end
  end
  return nothing
end

function init_vector!(v::Array{T, 1}, mesh::CubicMeshGPU, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx, j = 1:mesh.ny, k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        vec_value = init(i,j,k, mesh.nx, mesh.ny, mesh.nz)
        if vec_value != nothing
          b[:, id] .= vec_value[:]
        end
  end
end

function init_vector!(v::Array{T, 1}, mesh::CylindricalTubeMeshGPU, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nxyz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nr, k = 1:mesh.nz
    id = index(i, 1, k, mesh.nr, 1, mesh.nz)
    b[:, id] .=  init(i, k, mesh.nr, mesh.nz)
  end
end

function init_vector!(v::CuArray{T, 1}, mesh::Mesh, init::Any) where {T<:AbstractFloat}
    F = _cuda_using_double.x ? Float64 : Float32
    tmp = zeros(F, 3*mesh.nxyz)
    init_vector!(tmp, mesh, init)
    copyto!(v, tmp)
    return true
end

function omega_to_spin(omega::CuArray{T, 1}, spin::CuArray{T, 1}, spin_next::CuArray{T, 1}, N::Int64) where {T<:AbstractFloat}
  #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
  #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
  function __kernal!(a, b, c, n::Int64)
      i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
      if 0 < i <= n
        j = 3*i-2
        @inbounds w1 = a[j]*0.5
        @inbounds w2 = a[j+1]*0.5
        @inbounds w3 = a[j+2]*0.5
        @inbounds m1 = b[j]
        @inbounds m2 = b[j+1]
        @inbounds m3 = b[j+2]
        r = 1 + w1*w1 + w2*w2 + w3*w3
        a11 = 1 + w1*w1 - w2*w2 - w3*w3
        a12 = 2*(w1*w2 - w3)
        a13 = 2*(w2 + w1*w3)
        a21 = 2*(w1*w2 + w3)
        a22 = 1 - w1*w1 + w2*w2 - w3*w3
        a23 = -2*(w1-w2*w3)
        a31 = 2*(-w2+w1*w3)
        a32 = 2*(w1+w2*w3)
        a33 = 1 - w1*w1 - w2*w2 + w3*w3
        @inbounds c[j] = (a11*m1 + a12*m2 + a13*m3)/r
        @inbounds c[j+1] = (a21*m1 + a22*m2 + a23*m3)/r
        @inbounds c[j+2] = (a31*m1 + a32*m2 + a33*m3)/r
      end
      return nothing
  end
  blk, thr = cudims(N)
  @cuda blocks=blk threads=thr __kernal!(omega, spin, spin_next, N)
  return nothing
end

function compute_dm!(dm::CuArray{T, 1}, m1::CuArray{T, 1}, m2::CuArray{T, 1}, N::Int64) where {T<:AbstractFloat}
  function __kernal!(c, a, b, n::Int64)
     i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
     if 0 < i <= n
         j = 3*i - 2
         @inbounds mx = a[j] - b[j]
         @inbounds my = a[j+1] - b[j+1]
         @inbounds mz = a[j+2] - b[j+2]
         @inbounds c[i] = CUDA.sqrt(mx*mx + my*my + mz*mz)
     end
     return nothing
  end
  blk, thr = cudims(N)
  @cuda blocks=blk threads=thr __kernal!(dm, m1, m2, N)
  return nothing
end


function normalise(a::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}
    function __kernal!(a,  n::Int64)
       i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
       if 0 < i <= n
           j = 3*i - 2
		   @inbounds m2 = a[j]*a[j] + a[j+1]*a[j+1] + a[j+2]*a[j+2]
           if m2>0
               @inbounds length = CUDA.rsqrt(m2)
               @inbounds a[j] *= length
               @inbounds a[j+1] *= length
               @inbounds a[j+2] *= length
           end
       end
       return nothing
    end
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr __kernal!(a, N)

   return nothing
end


function average_m(sim::AbstractSimGPU)
  b = reshape(sim.spin, 3, sim.nxyz)
  return Tuple(sum(b, dims=2)./sim.nxyz)  #FIXME: using non-zero N in general case?
end

function average_m(sim::AbstractSimGPU, shape::CuArray{Bool, 1})
  b = reshape(sim.spin, 3, sim.nxyz)
  return Tuple(sum(b, dims=2)./sum(shape))
end

function average_m_with_shape(sim::AbstractSimGPU)
    shape = sim.shape
  b = reshape(sim.spin, 3, sim.nxyz)
  return Tuple(sum(b, dims=2)./sum(shape))
end

function compute_skyrmion_number(v::Array{T, 1}, m::Array{T, 1}, mesh::TriangularMeshGPU) where {T<:AbstractFloat}
    nx,ny = mesh.nx, mesh.ny
    ngbs = Array(mesh.ngbs)
    for j = 1:ny, i=1:nx
        id = index(i, j, nx, ny)
        v[id] = T(0)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)

        id1 = 3*ngbs[1, id]
        id2 = 3*ngbs[2, id]
        if id1>0 && id2>0
                sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
                sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
                v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*ngbs[4, id]
        id2 = 3*ngbs[5, id]
        if id1>0 && id2>0
                sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
                sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
                v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4*pi);
    end
    return nothing
end


function compute_skyrmion_number(m::Array{T, 1}, mesh::TriangularMeshGPU) where {T<:AbstractFloat}
    nx,ny = mesh.nx, mesh.ny
    v = zeros(T, nx*ny)
    compute_skyrmion_number(v, m, mesh)
    return sum(v)
end


function compute_skyrmion_number(m::Array{T, 1}, mesh::CylindricalTubeMeshGPU) where {T<:AbstractFloat}
    dx = 2*pi*mesh.R/mesh.nr
    mm = convert_m_to_cylindrical(m, mesh.nr, mesh.nz)
    v = zeros(T, mesh.nr*mesh.nz)
    pbc_z = mesh.zperiodic ? "z" : ""
    mesh = FDMesh(nx=mesh.nr, ny=mesh.nz, nz=1, dx = dx, dy=mesh.dz, dz=1e-9, pbc="x"*pbc_z)
    compute_skyrmion_number(v, mm, mesh)
    return sum(v)
end

"""
    compute_guiding_center(sim::AbstractSim, mesh::CylindricalTubeMeshGPU)

compute the guiding center when a CylindricalTubeMeshGPU is used. However, the direct calculation will give 
an incorrect results if the skyrmion touches the edges. Since periodic boundary conditions are naturally 
imposed on the mesh, we can shift the spin forward and backward one third and recompute the guiding centers.
In principle, two of those three centers should always be the same no matter where the skyrmion is. 
"""
function compute_guiding_center(sim::AbstractSim, mesh::CylindricalTubeMeshGPU) where {T<:AbstractFloat}
  
  # two of the three numbers are very close, so we compute the distance between them, 
  # and return the one with the least distance
  function find_value_min_distance(x, y, z)
     all = [x,y,z]
     v = argmin(abs.(circshift(all, 1) - all))
     return all[v]
  end

    spin = Array(sim.spin)
    dx = 2*pi*mesh.R/mesh.nr
    dy = mesh.dz # we set dy using mesh.dz
    m = convert_m_to_cylindrical(spin, mesh.nr, mesh.nz) # it seems we don't have to transform it?
    mesh = FDMesh(nx=mesh.nr, ny=mesh.nz, nz=1, dx = dx, dy = dy, dz=1e-9, pbc="xy")
    Rx, Ry = compute_guiding_center(m, mesh, z=1)

    roll_nx = Int(div(mesh.nx, 3)) 
    roll_ny = Int(div(mesh.ny, 3))

    m2 = reshape(m, 3, mesh.nx, mesh.ny)
    m2f = circshift(m2, (0, roll_nx, roll_ny)) # roll m and recompute guiding center
    m2f = reshape(m2f, 3*mesh.nx*mesh.ny)
    Rxf, Ryf = compute_guiding_center(m2f, mesh, z=1)
    Rxf = (Rxf - roll_nx*mesh.dx < 0) ? (Rxf + (mesh.nx - roll_nx)*mesh.dx) : Rxf - roll_nx*mesh.dx
    Ryf = (Ryf - roll_ny*mesh.dy < 0) ? (Ryf + (mesh.ny - roll_ny)*mesh.dy) : Ryf - roll_ny*mesh.dy
    
    m2b = circshift(m2, (0, -roll_nx, -roll_ny)) 
    m2b = reshape(m2b, 3*mesh.nx*mesh.ny)
    Rxb, Ryb = compute_guiding_center(m2b, mesh, z=1)
    Rxb = (Rxb + roll_nx*mesh.dx)%(mesh.nx*dx)
    Ryb = (Ryb + roll_ny*mesh.dy)%(mesh.ny*dy)
    
    Rx = find_value_min_distance(Rx, Rxf, Rxb)
    Ry = find_value_min_distance(Ry, Ryf, Ryb)
    return Rx, Ry

end


