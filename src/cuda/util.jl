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
  for i = 1:mesh.nx
    for j = 1:mesh.ny
        for k = 1:mesh.nz
            id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
            b[:, id] .=  init(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        end
      end
  end
  return nothing
end

function init_vector!(v::Array{T, 1}, mesh::CubicMeshGPU, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
      for k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        b[:, id] .=  init(i,j,k, mesh.nx, mesh.ny, mesh.nz)
      end
    end
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
