using CUDAnative, CuArrays

function zeeman_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                        h_static::CuDeviceArray{T, 1}, energy::CuDeviceArray{T, 1},
                        Ms::CuDeviceArray{T, 1}, volume::T, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        mu0 = 4*pi*1e-7
        j = 3*(index-1)
        @inbounds h[j+1] = h_static[j+1]
        @inbounds h[j+2] = h_static[j+2]
        @inbounds h[j+3] = h_static[j+3]
        @inbounds energy[index] = -mu0*Ms[index]*volume*(m[j+1]*h[j+1] + m[j+2]*h[j+2] + m[j+3]*h[j+3])
    end
   return nothing
end


function time_zeeman_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                        h_static::CuDeviceArray{T, 1}, energy::CuDeviceArray{T, 1},
                        Ms::CuDeviceArray{T, 1}, volume::T, fx::T, fy::T, fz::T, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        mu0 = 4*pi*1e-7
        j = 3*(index-1)
        @inbounds h[j+1] = h_static[j+1]*fx
        @inbounds h[j+2] = h_static[j+2]*fy
        @inbounds h[j+3] = h_static[j+3]*fz
        @inbounds energy[index] = -mu0*Ms[index]*volume*(m[j+1]*h[j+1] + m[j+2]*h[j+2] + m[j+3]*h[j+3])
    end
   return nothing
end

function exchange_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                     energy::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1}, A::T,
                     dx::T, dy::T, dz::T, nx::Int64, ny::Int64, nz::Int64, volume::T,
                     xperiodic::Bool, yperiodic::Bool, zperiodic::Bool) where {T<:AbstractFloat}

  index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

  nxy = nx * ny
  nxyz = nxy * nz
  ax = 2 * A / (dx * dx)
  ay = 2 * A / (dy * dy)
  az = 2 * A / (dz * dz)

  if 0 < index <= nxyz
      i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])

      indexm = 3*index
      @inbounds mx = m[indexm-2]
      @inbounds my = m[indexm-1]
      @inbounds mz = m[indexm]
      @inbounds Ms_local = Ms[index]

      if Ms_local == T(0)
          @inbounds energy[index] = 0
          @inbounds h[indexm - 2] = 0
          @inbounds h[indexm - 1] = 0
          @inbounds h[indexm] = 0
          return nothing
      end
      id, idm = 0, 0
      fx, fy, fz = T(0.0), T(0.0), T(0.0)
      if k>1 || zperiodic
          id = (k==1) ? index - nxy + nxyz : index - nxy
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += az*(m[idm-2] - mx)
              @inbounds fy += az*(m[idm-1] - my)
              @inbounds fz += az*(m[idm] - mz)
          end
      end

      if j>1 || yperiodic
          id = (j==1) ? index - nx + nxy : index - nx
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += ay*(m[idm-2] - mx)
              @inbounds fy += ay*(m[idm-1] - my)
              @inbounds fz += ay*(m[idm] - mz)
         end
      end

      if i>1 || xperiodic
          id = (i==1) ? index - 1 + nx : index -1
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += ax*(m[idm-2] - mx)
              @inbounds fy += ax*(m[idm-1] - my)
              @inbounds fz += ax*(m[idm] - mz)
          end
      end

      if i<nx || xperiodic
          id = (i==nx) ? index +1 - nx : index +1
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += ax*(m[idm-2]- mx)
              @inbounds fy += ax*(m[idm-1]- my)
              @inbounds fz += ax*(m[idm]- mz)
          end
      end

      if j<ny || yperiodic
          id = (j==ny) ? index + nx - nxy : index + nx
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += ay*(m[idm-2]- mx)
              @inbounds fy += ay*(m[idm-1]- my)
              @inbounds fz += ay*(m[idm]- mz)
          end
      end

      if k<nz || zperiodic
          id = (k==nz) ? index + nxy - nxyz : index + nxy
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += az*(m[idm-2]- mx)
              @inbounds fy += az*(m[idm-1]- my)
              @inbounds fz += az*(m[idm]- mz)
          end
      end
      Ms_inv = 1.0/(4.0*pi*1e-7*Ms_local)
      @inbounds energy[index] = - 0.5 * (fx * mx + fy * my + fz * mz)* volume;
      @inbounds h[indexm - 2] =  fx*Ms_inv;
      @inbounds h[indexm - 1] =  fy*Ms_inv;
      @inbounds h[indexm] = fz*Ms_inv;
  end
  return nothing
end


function bulkdmi_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                     energy::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1},
                     Dx::T, Dy::T, Dz::T, dx::T, dy::T, dz::T, nx::Int64,
                     ny::Int64, nz::Int64, volume::T, xperiodic::Bool,
                     yperiodic::Bool, zperiodic::Bool) where {T<:AbstractFloat}

  index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

  nxy = nx * ny
  nxyz = nxy * nz

  if 0< index <= nxyz
      i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])

      indexm = 3*index
      @inbounds mx = m[indexm-2]
      @inbounds my = m[indexm-1]
      @inbounds mz = m[indexm]
      @inbounds Ms_local = Ms[index]

      if Ms_local == T(0)
          @inbounds energy[index] = 0
          @inbounds h[indexm - 2] = 0
          @inbounds h[indexm - 1] = 0
          @inbounds h[indexm] = 0
          return nothing
      end
      id, idm = 0, 0
      fx, fy, fz = T(0.0), T(0.0), T(0.0)
      if k>1 || zperiodic
          id = (k==1) ? index - nxy + nxyz : index - nxy
          idm = 3*id - 2
		  @inbounds fx += (Dz/dz)*cross_x(T(0),T(0),T(1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (Dz/dz)*cross_y(T(0),T(0),T(1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (Dz/dz)*cross_z(T(0),T(0),T(1),m[idm],m[idm+1],m[idm+2]);
      end

      if j>1 || yperiodic
          id = (j==1) ? index - nx + nxy : index - nx
		  idm = 3*id - 2
		  @inbounds fx += (Dy/dy)*cross_x(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (Dy/dy)*cross_y(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (Dy/dy)*cross_z(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if i>1 || xperiodic
          id = (i==1) ? index - 1 + nx : index -1
          idm = 3*id - 2
		  @inbounds fx += (Dx/dx)*cross_x(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (Dx/dx)*cross_y(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (Dx/dx)*cross_z(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if i<nx || xperiodic
          id = (i==nx) ? index +1 - nx : index +1
          idm = 3*id - 2
		  @inbounds fx += (Dx/dx)*cross_x(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (Dx/dx)*cross_y(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (Dx/dx)*cross_z(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if j<ny || yperiodic
          id = (j==ny) ? index + nx - nxy : index + nx
          idm = 3*id - 2
		  @inbounds fx += (Dy/dy)*cross_x(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (Dy/dy)*cross_y(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (Dy/dy)*cross_z(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if k<nz || zperiodic
          id = (k==nz) ? index + nxy - nxyz : index + nxy
		  idm = 3*id - 2
		  @inbounds fx += (Dz/dz)*cross_x(T(0),T(0),T(-1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (Dz/dz)*cross_y(T(0),T(0),T(-1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (Dz/dz)*cross_z(T(0),T(0),T(-1),m[idm],m[idm+1],m[idm+2]);
      end
      Ms_inv = 1.0/(4.0*pi*1e-7*Ms_local)
      @inbounds energy[index] = - 0.5 * (fx * mx + fy * my + fz * mz)* volume;
      @inbounds h[indexm - 2] =  fx*Ms_inv;
      @inbounds h[indexm - 1] =  fy*Ms_inv;
      @inbounds h[indexm] = fz*Ms_inv;
  end
  return nothing
end

function interfacial_dmi_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                     energy::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1},
                     D::T, dx::T, dy::T, dz::T, nx::Int64,
                     ny::Int64, nz::Int64, volume::T, xperiodic::Bool,
                     yperiodic::Bool) where {T<:AbstractFloat}

  index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

  nxy = nx * ny
  nxyz = nxy * nz

  if 0< index <= nxyz
      i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])

      indexm = 3*index
      @inbounds mx = m[indexm-2]
      @inbounds my = m[indexm-1]
      @inbounds mz = m[indexm]
      @inbounds Ms_local = Ms[index]

      if Ms_local == T(0)
          @inbounds energy[index] = 0
          @inbounds h[indexm - 2] = 0
          @inbounds h[indexm - 1] = 0
          @inbounds h[indexm] = 0
          return nothing
      end
      id, idm = 0, 0
      fx, fy, fz = T(0.0), T(0.0), T(0.0)

      if j>1 || yperiodic
          id = (j==1) ? index - nx + nxy : index - nx
		  idm = 3*id - 2
		  @inbounds fx += (D/dy)*cross_x(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dy)*cross_y(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dy)*cross_z(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if i>1 || xperiodic
          id = (i==1) ? index - 1 + nx : index -1
          idm = 3*id - 2
		  @inbounds fx += (D/dx)*cross_x(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dx)*cross_y(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dx)*cross_z(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if i<nx || xperiodic
          id = (i==nx) ? index +1 - nx : index +1
          idm = 3*id - 2
		  @inbounds fx += (D/dx)*cross_x(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dx)*cross_y(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dx)*cross_z(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if j<ny || yperiodic
          id = (j==ny) ? index + nx - nxy : index + nx
          idm = 3*id - 2
		  @inbounds fx += (D/dy)*cross_x(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dy)*cross_y(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dy)*cross_z(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      Ms_inv = 1.0/(4.0*pi*1e-7*Ms_local)
      @inbounds energy[index] = - 0.5 * (fx * mx + fy * my + fz * mz)* volume;
      @inbounds h[indexm - 2] =  fx*Ms_inv;
      @inbounds h[indexm - 1] =  fy*Ms_inv;
      @inbounds h[indexm] = fz*Ms_inv;
  end
  return nothing
end


function spatial_bulkdmi_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                     energy::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1},
                     Ds::CuDeviceArray{T, 1},dx::T, dy::T, dz::T, nx::Int64,
                     ny::Int64, nz::Int64, volume::T, xperiodic::Bool,
                     yperiodic::Bool, zperiodic::Bool) where {T<:AbstractFloat}

  index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

  nxy = nx * ny
  nxyz = nxy * nz

  if 0< index <= nxyz
      i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])

      indexm = 3*index
      @inbounds mx = m[indexm-2]
      @inbounds my = m[indexm-1]
      @inbounds mz = m[indexm]
      @inbounds Ms_local = Ms[index]
      @inbounds D = Ds[index]

      if Ms_local == T(0)
          @inbounds energy[index] = 0
          @inbounds h[indexm - 2] = 0
          @inbounds h[indexm - 1] = 0
          @inbounds h[indexm] = 0
          return nothing
      end
      id, idm = 0, 0
      fx, fy, fz = T(0.0), T(0.0), T(0.0)
      if k>1 || zperiodic
          id = (k==1) ? index - nxy + nxyz : index - nxy
          idm = 3*id - 2
		  @inbounds fx += (D/dz)*cross_x(T(0),T(0),T(1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dz)*cross_y(T(0),T(0),T(1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dz)*cross_z(T(0),T(0),T(1),m[idm],m[idm+1],m[idm+2]);
      end

      if j>1 || yperiodic
          id = (j==1) ? index - nx + nxy : index - nx
		  idm = 3*id - 2
		  @inbounds fx += (D/dy)*cross_x(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dy)*cross_y(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dy)*cross_z(T(0),T(1),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if i>1 || xperiodic
          id = (i==1) ? index - 1 + nx : index -1
          idm = 3*id - 2
		  @inbounds fx += (D/dx)*cross_x(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dx)*cross_y(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dx)*cross_z(T(1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if i<nx || xperiodic
          id = (i==nx) ? index +1 - nx : index +1
          idm = 3*id - 2
		  @inbounds fx += (D/dx)*cross_x(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dx)*cross_y(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dx)*cross_z(T(-1),T(0),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if j<ny || yperiodic
          id = (j==ny) ? index + nx - nxy : index + nx
          idm = 3*id - 2
		  @inbounds fx += (D/dy)*cross_x(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dy)*cross_y(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dy)*cross_z(T(0),T(-1),T(0),m[idm],m[idm+1],m[idm+2]);
      end

      if k<nz || zperiodic
          id = (k==nz) ? index + nxy - nxyz : index + nxy
		  idm = 3*id - 2
		  @inbounds fx += (D/dz)*cross_x(T(0),T(0),T(-1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fy += (D/dz)*cross_y(T(0),T(0),T(-1),m[idm],m[idm+1],m[idm+2]);
		  @inbounds fz += (D/dz)*cross_z(T(0),T(0),T(-1),m[idm],m[idm+1],m[idm+2]);
      end
      Ms_inv = 1.0/(4.0*pi*1e-7*Ms_local)
      @inbounds energy[index] = - 0.5 * (fx * mx + fy * my + fz * mz)* volume;
      @inbounds h[indexm - 2] =  fx*Ms_inv;
      @inbounds h[indexm - 1] =  fy*Ms_inv;
      @inbounds h[indexm] = fz*Ms_inv;
  end
  return nothing
end

function anisotropy_kernel!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                        energy::CuDeviceArray{T, 1}, Ku::CuDeviceArray{T, 1},
                        axis_x::T, axis_y::T, axis_z::T,
                        Ms::CuDeviceArray{T, 1}, volume::T, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        mu0 = 4*pi*1e-7
        j = 3*(index-1)
        @inbounds Ms_local = Ms[index]
        if Ms_local == 0.0
            return nothing
        end
        Ms_inv::T = 2.0/(mu0*Ms_local)
        @inbounds sa = m[j+1]*axis_x+m[j+2]*axis_y+m[j+3]*axis_z
        @inbounds h[j+1] = Ku[index]*m[j+1]*axis_x*Ms_inv
        @inbounds h[j+2] = Ku[index]*m[j+2]*axis_y*Ms_inv
        @inbounds h[j+3] = Ku[index]*m[j+3]*axis_z*Ms_inv
        @inbounds energy[index] = Ku[index]*(1.0-sa*sa)*volume
    end
   return nothing
end


function exchange_kernel_rkky!(m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                        energy::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1},
                        sigma::T, nx::Int64, ny::Int64, nz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nx*ny*nz
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
        if k!=1
            return nothing
        end

        id1 = (j-1) * nx + i
        id2 = (nz-1) * nx*ny + (j-1) * nx + i

        k1 = 3*id1-2
        k2 = 3*id2-2
        @inbounds mbx = m[k1]
        @inbounds mby = m[k1+1]
        @inbounds mbz = m[k1+2]

        @inbounds mtx = m[k2]
        @inbounds mty = m[k2+1]
        @inbounds mtz = m[k2+2]

        mu0 = 4*pi*1e-7

        @inbounds Ms1 =  Ms[id1]
        @inbounds Ms2 =  Ms[id2]
        if Ms1> 0 && Ms2 > 0
            Ms_inv = 1.0/(Ms1*mu0)
            @inbounds h[k1] = sigma*Ms_inv*mtx
            @inbounds h[k1+1] = sigma*Ms_inv*mty
            @inbounds h[k1+2] = sigma*Ms_inv*mtz
            @inbounds energy[id1] = -0.5*sigma*(1-mtx*mbx-mty*mby-mtz*mbz)

            Ms_inv = 1.0/(Ms2*mu0)
            @inbounds h[k2] = sigma*Ms_inv*mbx
            @inbounds h[k2+1] = sigma*Ms_inv*mby
            @inbounds h[k2+2] = sigma*Ms_inv*mbz
            @inbounds energy[id2] = energy[id1]
        end

    end
   return nothing
end
