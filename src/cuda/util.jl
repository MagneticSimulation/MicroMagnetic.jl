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


function  init_scalar!(v::CuArray{T, 1}, mesh::Mesh, init::Number) where {T<:AbstractFloat}
	v[:] .= init
    return nothing
end

function init_scalar!(v::CuArray{T, 1}, mesh::Mesh, init_fun::Function) where {T<:AbstractFloat}
    mesh = sim.mesh
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        v[id] = init_fun(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    return true
end

function omega_to_spin(omega::CuArray{T, 1}, spin::CuArray{T, 1}, spin_next::CuArray{T, 1}, N::Int64) where {T<:AbstractFloat}
  #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
  #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
  function __kernal!(a, b, c, N::Int64)
      i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
      if 0 < i <= N
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
  blk, thr = CuArrays.cudims(N)
  @cuda blocks=blk threads=thr __kernal!(omega, spin, spin_next, N)
  return nothing
end


function compute_dmdt(m1::CuArray{T, 1}, m2::CuArray{T, 1}, N::Int64, dt::Float64) where {T<:AbstractFloat}
  dmdt = cuzeros(T, N)
  function __kernal!(a, b, c, n::Int64)
     i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
     if 0 < i <= n
         j = 3*i - 2
         @inbounds mx = a[j] - b[j]
         @inbounds my = a[j+1] - b[j+1]
         @inbounds mz = a[j+2] - b[j+2]
         @inbounds c[i] = CUDAnative.sqrt(mx*mx + my*my + mz*mz)
     end
     return nothing
  end
  blk, thr = CuArrays.cudims(N)
  @cuda blocks=blk threads=thr __kernal!(m1, m2, dmdt, N)
  max_dmdt = maximum(dmdt)/dt #TODO: (1)how to free dmdt? (2)how to avoid to use maximum?
  return max_dmdt
end


function normalise(a::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}
    function __kernal!(a,  n::Int64)
       i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
       if 0 < i <= n
           j = 3*i - 2
           @inbounds length = 1.0/CUDAnative.sqrt(a[j]^2 + a[j+1]^2 + a[j+2]^2)
           @inbounds a[j] *= length
           @inbounds a[j+1] *= length
           @inbounds a[j+2] *= length
       end
       return nothing
    end
    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr __kernal!(a, N)

   return nothing
end


function average_m(sim::AbstractSimGPU)
  b = reshape(sim.spin, 3, sim.nxyz)
  return Tuple(sum(b, dims=2)./sim.nxyz)  #FIXME: using non-zero N in general case?
end
