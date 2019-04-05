function init_vector!(v::Array{Float64, 1}, mesh::FDMeshGPU, init::Function)
  nxyz = mesh.nx*mesh.ny*mesh.nz
  dx,dy,dz = mesh.dx, mesh.dy, mesh.dz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
      for k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        b[:, id] .=  init(i,j,k,dx,dy,dz)
      end
    end
  end
  return nothing
end

function omega_to_spin(omega::CuArray{T, 1}, spin::CuArray{T, 1}, spin_next::CuArray{T, 1}, N::Int64) where {T<:AbstractFloat}
  #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
  #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
  function kernal!(a, b, c, N::Int64)
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
  end
  blk, thr = CuArrays.cudims(N)
  @cuda blocks=blk threads=thr kernel!(omega, spin, spin_next, N)
  return nothing
end
