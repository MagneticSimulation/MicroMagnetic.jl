using CUDAdrv, CUDAnative, CuArrays

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
        @inbounds energy[i] = -mu0*Ms[i]*volume*(spin[j+1]*h[j+1] + spin[j+2]*h[j+2] + spin[j+3]*h[j+3])
    end
   return nothing
end

function effective_field(zeeman::ZeemanGPU, sim::MicroSimGPU, spin::CuArray{FloatGPU, 1}, t::Float64)
  mu0 = 4*pi*1e-7
  nxyz = sim.nxyz
  volume = sim.mesh.volume
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n zeeman_kernel!(spin, sim.field, zeeman.field_gpu, energy, sim.Ms, volume, nxyz)
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

  if 0< index <= nxyz
      i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])

      indexm = 3*index
      @inbounds mx = m[indexm-2]
      @inbounds my = m[indexm-1]
      @inbounds mz = m[indexm]
      @inbounds Ms_local = Ms[index]

      if Ms_local == 0.0
          return nothing
      end
      id, idm = 0, 0
      fx, fy, fz = FloatGPU(0.0), FloatGPU(0.0), FloatGPU(0.0)
      if k>1 || zperiodic
          id = (k==1) ? index - nxy + nxyz : index - nxy
          idm = 3*id
          @inbounds fx += az*(m[idm-2]- mx)
          @inbounds fy += az*(m[idm-1]- my)
          @inbounds fz += az*(m[idm]- mz)
      end

      if j>1 || yperiodic
          id = (j==1) ? index - nx + nxy : index - nx
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += ay*(m[idm-2]- mx)
              @inbounds fy += ay*(m[idm-1]- my)
              @inbounds fz += ay*(m[idm]- mz)
         end
      end

      if i>1 || xperiodic
          id = (i==1) ? index - 1 + nx : index -1
          idm = 3*id
          @inbounds if Ms[id]>0
              @inbounds fx += ax*(m[idm-2]- mx)
              @inbounds fy += ax*(m[idm-1]- my)
              @inbounds fz += ax*(m[idm]- mz)
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

function effective_field(exch::ExchangeGPU, sim::MicroSimGPU, spin::CuArray{FloatGPU, 1}, t::Float64)
  mu0 = 4*pi*1e-7
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_kernel!(spin, sim.field, sim.energy, sim.Ms, exch.A,
                    mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz, mesh.volume,
                    mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end


function effective_field(sim::MicroSimGPU, spin::CuArray{FloatGPU, 1}, t::Float64)
  fill!(sim.driver.field, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.driver.field .+= sim.field
  end
  return 0
end
