using CuArrays.CUFFT
using FFTW
using LinearAlgebra

function newell_f_gpu(x::Float64, y::Float64, z::Float64)::Float64
   x2 = x*x;
   y2 = y*y;
   z2 = z*z;
   R = CUDAnative.sqrt(x2+y2+z2);
   if R == 0.0
     return 0.0;
   end

   f = 1.0/6*(2*x2-y2-z2)*R

   if x2>0
    f -= x*y*z*CUDAnative.atan(y*z/(x*R));
   end

   if x2+z2>0
     f += 0.5*y*(z2-x2)*CUDAnative.asinh(y/(CUDAnative.sqrt(x2+z2)))
   end

   if x2+y2>0
     f += 0.5*z*(y2-x2)*CUDAnative.asinh(z/(CUDAnative.sqrt(x2+y2)))
   end
  return f;
end


function newell_g_gpu(x::Float64, y::Float64, z::Float64)::Float64
   x2 = x*x;
   y2 = y*y;
   z2 = z*z;

   R = CUDAnative.sqrt(x2+y2+z2);
   if R == 0.0
     return 0.0;
   end

   g = -1.0/3*x*y*R;

   if z2>0   g -= 1.0/6*z2*z*CUDAnative.atan(x*y/(z*R)) end
   if y2>0   g -= 0.5*y2*z*CUDAnative.atan(x*z/(y*R)) end
   if x2>0   g -= 0.5*x2*z*CUDAnative.atan(y*z/(x*R)) end

   if x2+y2>0
     g += x*y*z*CUDAnative.asinh(z/(CUDAnative.sqrt(x2+y2)))
    end

   if y2+z2>0
     g += 1.0/6*y*(3*z2-y2)*CUDAnative.asinh(x/(CUDAnative.sqrt(y2+z2)))
   end

   if x2+z2>0
     g += 1.0/6*x*(3*z2-x2)*CUDAnative.asinh(y/(CUDAnative.sqrt(x2+z2)))
   end

  return g;
end


function demag_tensor_xx_gpu(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)

  tensor = 8.0*newell_f_gpu(x,y,z);

  tensor -= 4.0*newell_f_gpu(x+dx,y,z);
  tensor -= 4.0*newell_f_gpu(x-dx,y,z);
  tensor -= 4.0*newell_f_gpu(x,y-dy,z);
  tensor -= 4.0*newell_f_gpu(x,y+dy,z);
  tensor -= 4.0*newell_f_gpu(x,y,z-dz);
  tensor -= 4.0*newell_f_gpu(x,y,z+dz);

  tensor +=  2.0*newell_f_gpu(x+dx,y+dy,z);
  tensor +=  2.0*newell_f_gpu(x+dx,y-dy,z);
  tensor +=  2.0*newell_f_gpu(x-dx,y-dy,z);
  tensor +=  2.0*newell_f_gpu(x-dx,y+dy,z);
  tensor +=  2.0*newell_f_gpu(x+dx,y,z+dz);
  tensor +=  2.0*newell_f_gpu(x+dx,y,z-dz);
  tensor +=  2.0*newell_f_gpu(x-dx,y,z+dz);
  tensor +=  2.0*newell_f_gpu(x-dx,y,z-dz);
  tensor +=  2.0*newell_f_gpu(x,y-dy,z-dz);
  tensor +=  2.0*newell_f_gpu(x,y-dy,z+dz);
  tensor +=  2.0*newell_f_gpu(x,y+dy,z+dz);
  tensor +=  2.0*newell_f_gpu(x,y+dy,z-dz);

  tensor -= newell_f_gpu(x+dx,y+dy,z+dz);
  tensor -= newell_f_gpu(x+dx,y+dy,z-dz);
  tensor -= newell_f_gpu(x+dx,y-dy,z+dz);
  tensor -= newell_f_gpu(x+dx,y-dy,z-dz);
  tensor -= newell_f_gpu(x-dx,y+dy,z+dz);
  tensor -= newell_f_gpu(x-dx,y+dy,z-dz);
  tensor -= newell_f_gpu(x-dx,y-dy,z+dz);
  tensor -= newell_f_gpu(x-dx,y-dy,z-dz);

  return tensor/(4.0*pi*dx*dy*dz)

end


function demag_tensor_xy_gpu(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)

  tensor = 8.0*newell_g_gpu(x,y,z);

  tensor -= 4.0*newell_g_gpu(x+dx,y,z);
  tensor -= 4.0*newell_g_gpu(x-dx,y,z);
  tensor -= 4.0*newell_g_gpu(x,y-dy,z);
  tensor -= 4.0*newell_g_gpu(x,y+dy,z);
  tensor -= 4.0*newell_g_gpu(x,y,z-dz);
  tensor -= 4.0*newell_g_gpu(x,y,z+dz);


  tensor +=  2.0*newell_g_gpu(x+dx,y+dy,z);
  tensor +=  2.0*newell_g_gpu(x+dx,y-dy,z);
  tensor +=  2.0*newell_g_gpu(x-dx,y-dy,z);
  tensor +=  2.0*newell_g_gpu(x-dx,y+dy,z);
  tensor +=  2.0*newell_g_gpu(x+dx,y,z+dz);
  tensor +=  2.0*newell_g_gpu(x+dx,y,z-dz);
  tensor +=  2.0*newell_g_gpu(x-dx,y,z+dz);
  tensor +=  2.0*newell_g_gpu(x-dx,y,z-dz);
  tensor +=  2.0*newell_g_gpu(x,y-dy,z-dz);
  tensor +=  2.0*newell_g_gpu(x,y-dy,z+dz);
  tensor +=  2.0*newell_g_gpu(x,y+dy,z+dz);
  tensor +=  2.0*newell_g_gpu(x,y+dy,z-dz);

  tensor -= newell_g_gpu(x+dx,y+dy,z+dz);
  tensor -= newell_g_gpu(x+dx,y+dy,z-dz);
  tensor -= newell_g_gpu(x+dx,y-dy,z+dz);
  tensor -= newell_g_gpu(x+dx,y-dy,z-dz);
  tensor -= newell_g_gpu(x-dx,y+dy,z+dz);
  tensor -= newell_g_gpu(x-dx,y+dy,z-dz);
  tensor -= newell_g_gpu(x-dx,y-dy,z+dz);
  tensor -= newell_g_gpu(x-dx,y-dy,z-dz);

  return tensor/(4.0*pi*dx*dy*dz)
end

function compute_tensors_kernel_xx!(tensor_xx, tensor_yy, tensor_zz,
                                 nx::Int64, ny::Int64, nz::Int64,
                                 dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(tensor_xx)
    if 0 < index <= nx_fft*ny_fft*nz_fft
        i,j,k = Tuple(CuArrays.CartesianIndices(tensor_xx)[index])
        if (nx_fft%2 == 0 && i == nx+1) || (ny_fft%2 == 0 && j == ny+1) || (nz_fft%2 == 0 && k == nz+1)
          return nothing
        end
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
            x = (i<=nx) ? (i-1+p*nx)*dx : (i-nx_fft-1-p*nx)*dx
            y = (j<=ny) ? (j-1+q*ny)*dy : (j-ny_fft-1-q*ny)*dy
            z = (k<=nz) ? (k-1+s*nz)*dz : (k-nz_fft-1-s*nz)*dz
            @inbounds tensor_xx[i,j,k] = demag_tensor_xx_gpu(x,y,z,dx,dy,dz)
            @inbounds tensor_yy[i,j,k] = demag_tensor_xx_gpu(y,x,z,dy,dx,dz)
            @inbounds tensor_zz[i,j,k] = demag_tensor_xx_gpu(z,y,x,dz,dy,dx)
        end
    end
    return nothing
end

function compute_tensors_kernel_xy!(tensor_xy, tensor_xz, tensor_yz,
                                 nx::Int64, ny::Int64, nz::Int64,
                                 dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(tensor_xy)
    if 0 < index <= nx_fft*ny_fft*nz_fft
        i,j,k = Tuple(CuArrays.CartesianIndices(tensor_xy)[index])
        if (nx_fft%2 == 0 && i == nx+1) || (ny_fft%2 == 0 && j == ny+1) || (nz_fft%2 == 0 && k == nz+1)
          return nothing
        end
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
            x = (i<=nx) ? (i-1+p*nx)*dx : (i-nx_fft-1-p*nx)*dx
            y = (j<=ny) ? (j-1+q*ny)*dy : (j-ny_fft-1-q*ny)*dy
            z = (k<=nz) ? (k-1+s*nz)*dz : (k-nz_fft-1-s*nz)*dz
            @inbounds tensor_xy[i,j,k] = demag_tensor_xy_gpu(x,y,z,dx,dy,dz)
            @inbounds tensor_xz[i,j,k] = demag_tensor_xy_gpu(x,z,y,dx,dz,dy)
            @inbounds tensor_yz[i,j,k] = demag_tensor_xy_gpu(y,z,x,dy,dz,dx)
        end
    end
    return nothing
end

function distribute_m_kernel!(m_gpu, mx_gpu, my_gpu, mz_gpu, Ms, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(mx_gpu)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
		p = 3*index-2
        @inbounds mx_gpu[i,j,k] = m_gpu[p]*Ms[index]
        @inbounds my_gpu[i,j,k] = m_gpu[p+1]*Ms[index]
        @inbounds mz_gpu[i,j,k] = m_gpu[p+2]*Ms[index]
    end
    return nothing
end

function distribute_m_kernel_II!(m_gpu, mx_gpu, my_gpu, mz_gpu, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(mx_gpu)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
        mx_gpu[i,j,k] = m_gpu[3*index-2]
        my_gpu[i,j,k] = m_gpu[3*index-1]
        mz_gpu[i,j,k] = m_gpu[3*index]
    end
    return nothing
end

function collect_h_kernel!(h, hx, hy, hz, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0< index <= nx*ny*nz
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
		p = 3*index-2
        @inbounds h[p] = -1.0*hx[i,j,k]
        @inbounds h[p+1] = -1.0*hy[i,j,k]
        @inbounds h[p+2] = -1.0*hz[i,j,k]
    end
    return nothing
end


function compute_energy_kernel!(energy::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                               Ms::CuDeviceArray{T, 1}, volume::T, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0< index <= nxyz
		mu0 = 4*pi*1e-7
        j = 3*index-2
        @inbounds energy[index] = -0.5*mu0*Ms[index]*volume*(m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2])
    end
    return nothing
end
