function uniform_random_sphere_kernel!(m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        j = 3*index - 2
        @inbounds phi = rnd[j]*2*pi;
        @inbounds ct = 2*rnd[j+1] - 1;
        st = CUDA.sqrt(1-ct*ct);
        @inbounds m[j] = st*CUDA.cos(phi);
        @inbounds m[j+1] = st*CUDA.sin(phi);
        @inbounds m[j+2] = ct;
    end
   return nothing
end

function uniform_random_circle_xy_kernel!(m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        j = 3*index - 2
        @inbounds phi = rnd[j]*2*pi;
        @inbounds m[j] = CUDA.cos(phi);
        @inbounds m[j+1] = CUDA.sin(phi);
        @inbounds m[j+2] = 0;
    end
   return nothing
end
function run_monte_carlo_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1},
                              dE::CuDeviceArray{T, 1}, temp::Float64,
                              H::CuDeviceArray{T, 1},
                              J::CuDeviceArray{T, 1},
                              K::CuDeviceArray{T, 1},
                              TotRdN::Int64,neig::CuDeviceArray{Int32, 2},n_neig::CuDeviceArray{Int32, 1},
                              nx::Int64, ny::Int64, nz::Int64,bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nx*ny*nz #&& shape[index]
      a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
      # sign = cubic ? 1 : -1
      if mod(a+2*b+c, 5) != bias
          return nothing
      end
      i = 3*index - 2

      @inbounds mx=m[i]
      @inbounds my=m[i+1]
      @inbounds mz=m[i+2]
      @inbounds nmx=next_m[i]
      @inbounds nmy=next_m[i+1]
      @inbounds nmz=next_m[i+2]
      @inbounds dmx = nmx - mx
      @inbounds dmy = nmy - my
      @inbounds dmz = nmz - mz       

      delta_E = -(dmx*H[1]+dmy*H[2]+dmz*H[3]) #zeeman

      u1x=T(1.0)
      u2x=T(-0.5)
      u3x=T(-0.5)
      u1y=T(0.0)
      u2y=T(sqrt(3.0)/2.0)
      u3y=T(-sqrt(3.0)/2.0)

      mu1s=CUDA.pow(mx*u1x+my*u1y,2) 
      mu2s=CUDA.pow(mx*u2x+my*u2y,2) 
      mu3s=CUDA.pow(mx*u3x+my*u3y,2) 
      delta_E += K[1]*16/9.0*(mu1s*mu2s+mu2s*mu3s+mu3s*mu1s)+K[2]*16*(mu1s*mu2s*mu3s)
      nu1s=CUDA.pow(nmx*u1x+nmy*u1y,2) 
      nu2s=CUDA.pow(nmx*u2x+nmy*u2y,2) 
      nu3s=CUDA.pow(nmx*u3x+nmy*u3y,2)
      delta_E -= K[1]*16/9.0*(nu1s*nu2s+nu2s*nu3s+nu3s*nu1s)+K[2]*16*(nu1s*nu2s*nu3s)
      for iirdn = 1:TotRdN
        # nns=n_neig[iirdn]+1
        # nne=n_neig[iirdn+1]
      for nn = n_neig[iirdn]+1:n_neig[iirdn+1]
        id=neig[nn,index]
          if id>0
              k = 3*id-2
              @inbounds sx = m[k]
              @inbounds sy = m[k+1]
              @inbounds sz = m[k+2]
              delta_E -= J[iirdn]*(dmx*sx + dmy*sy + dmz*sz) #exchange
          end
      end
      end

      @inbounds dE[index] = delta_E
      update = false

      if delta_E < 0
          update = true
      else
          @inbounds rnd_i = rnd[i+2]
          if rnd_i < CUDA.exp(-delta_E/temp)
              update = true
          end
      end

      if update == true
        @inbounds m[i] = next_m[i];
        @inbounds m[i+1] = next_m[i+1];
        @inbounds m[i+2] = next_m[i+2]
      end
    end
   return nothing
end

function getEnergyAndStd!(m::CuDeviceArray{T, 1},m2::CuDeviceArray{T, 1},
            energy::CuDeviceArray{T, 1},energy2::CuDeviceArray{T, 1},
            qdens::CuDeviceArray{T, 1},H::CuDeviceArray{T, 1},
            J::CuDeviceArray{T, 1},K::CuDeviceArray{T, 1},
            TotRdN::Int64,neig::CuDeviceArray{Int32, 2},n_neig::CuDeviceArray{Int32, 1},
            nx::Int64, ny::Int64, nz::Int64,Volu::Bool) where {T<:AbstractFloat}

# @inline function Berg_Omega(ux::T, uy::T, uz::T, vx::T, vy::T, vz::T, wx::T, wy::T, wz::T) where {T<:AbstractFloat}
#     b = volume(ux, uy, uz, vx, vy, vz, wx, wy, wz)
#     a = 1.0 + (ux*vx + uy*vy + uz*vz) + (ux*wx + uy*wy + uz*wz) + (vx*wx + vy*wy + vz*wz)
#     return 2*CUDA.atan(b/a)
# end
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nx*ny*nz #&& shape[index]
      a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
      i = 3*index - 2

      @inbounds mx=m[i]
      @inbounds my=m[i+1]
      @inbounds mz=m[i+2]

      @inbounds En = -(mx*H[1]+my*H[2]+mz*H[3]) #zeeman

      u1x=T(1.0)
      u2x=T(-0.5)
      u3x=T(-0.5)
      u1y=T(0.0)
      u2y=T(sqrt(3.0)/2.0)
      u3y=T(-sqrt(3.0)/2.0)

      mu1s=CUDA.pow(mx*u1x+my*u1y,2)
      mu2s=CUDA.pow(mx*u2x+my*u2y,2) 
      mu3s=CUDA.pow(mx*u3x+my*u3y,2) 
      @inbounds En += -K[1]*16/9.0*(mu1s*mu2s+mu2s*mu3s+mu3s*mu1s)-K[2]*16*(mu1s*mu2s*mu3s)

      #第一近邻
      xperiodic=true
      yperiodic=true
      zperiodic=true
      
      for iirdn = 1:TotRdN
        @inbounds nns=n_neig[iirdn]+1
        @inbounds nne=n_neig[iirdn+1]
      for nn =nns:nne
        @inbounds id=neig[nn,index]
          if id>0
              k = 3*id-2
              @inbounds sx = m[k]
              @inbounds sy = m[k+1]
              @inbounds sz = m[k+2]
              @inbounds En += -0.5*J[iirdn]*(mx*sx + my*sy + mz*sz) #exchange
          end
      end
      end

      @inbounds energy[index] = En
      @inbounds energy2[index] = En*En
      #m2
      @inbounds m2[i]  =mx*mx
      @inbounds m2[i+1]=my*my
      @inbounds m2[i+2]=mz*mz
      #topological Dens
      sx1,sy1,sz1 = T(0),T(0),T(0)
      sx2,sy2,sz2 = T(0),T(0),T(0)
      # qdens=T(0)
      n11=3*neig[1,index]
      n12=3*neig[2,index]
      n14=3*neig[4,index]
      n15=3*neig[5,index]
      if n11>0 && n12>0
        @inbounds sx1,sy1,sz1 = m[n11-2],m[n11-1],m[n11]
        @inbounds sx2,sy2,sz2 = m[n12-2],m[n12-1],m[n12]
        # qtemp = Berg_Omega_Gpu(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1) 
        if Volu
          qtemp = volume(mx, my, mz, sx1, sy1, sz1, sx2, sy2, sz2 ) 
        else
          qtemp = Berg_Omega_Gpu(mx, my, mz, sx1, sy1, sz1, sx2, sy2, sz2 ) 
        end
      end
      if n14>0 && n15>0
        @inbounds sx1,sy1,sz1 = m[n14-2],m[n14-1],m[n14]
        @inbounds sx2,sy2,sz2 = m[n15-2],m[n15-1],m[n15]
        # qtemp += Berg_Omega_Gpu(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)  
        if Volu
          qtemp += volume(mx, my, mz, sx1, sy1, sz1, sx2, sy2, sz2 ) 
        else
          qtemp += Berg_Omega_Gpu(mx, my, mz, sx1, sy1, sz1, sx2, sy2, sz2 )  
        end
      end

      if Volu
        @inbounds qdens[index]=qtemp*0.75/pi
      else
        @inbounds qdens[index]=qtemp*0.25/pi
      end

    end
   return nothing
end


