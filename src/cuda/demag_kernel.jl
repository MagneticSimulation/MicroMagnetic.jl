using CUDA.CUFFT
using FFTW
using LinearAlgebra

function newell_f_gpu(x::Float64, y::Float64, z::Float64)::Float64
   x2 = x*x;
   y2 = y*y;
   z2 = z*z;
   R = CUDA.sqrt(x2+y2+z2);
   if R == 0.0
     return 0.0;
   end

   f = 1.0/6*(2*x2-y2-z2)*R

   if x2>0
    f -= x*y*z*CUDA.atan(y*z/(x*R));
   end

   if x2+z2>0
     f += 0.5*y*(z2-x2)*CUDA.asinh(y/(CUDA.sqrt(x2+z2)))
   end

   if x2+y2>0
     f += 0.5*z*(y2-x2)*CUDA.asinh(z/(CUDA.sqrt(x2+y2)))
   end
  return f;
end


function newell_g_gpu(x::Float64, y::Float64, z::Float64)::Float64
   x2 = x*x;
   y2 = y*y;
   z2 = z*z;

   R = CUDA.sqrt(x2+y2+z2);
   if R == 0.0
     return 0.0;
   end

   g = -1.0/3*x*y*R;

   if z2>0   g -= 1.0/6*z2*z*CUDA.atan(x*y/(z*R)) end
   if y2>0   g -= 0.5*y2*z*CUDA.atan(x*z/(y*R)) end
   if x2>0   g -= 0.5*x2*z*CUDA.atan(y*z/(x*R)) end

   if x2+y2>0
     g += x*y*z*CUDA.asinh(z/(CUDA.sqrt(x2+y2)))
    end

   if y2+z2>0
     g += 1.0/6*y*(3*z2-y2)*CUDA.asinh(x/(CUDA.sqrt(y2+z2)))
   end

   if x2+z2>0
     g += 1.0/6*x*(3*z2-x2)*CUDA.asinh(y/(CUDA.sqrt(x2+z2)))
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

function demag_tensor_yy_gpu(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xx_gpu(y,x,z,dy,dx,dz);
end

function demag_tensor_zz_gpu(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xx_gpu(z,y,x,dz,dy,dx);
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

function demag_tensor_xz_gpu(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xy_gpu(x,z,y,dx,dz,dy);
end

function demag_tensor_yz_gpu(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xy_gpu(y,z,x,dy,dz,dx)
end

function compute_tensors_kernel_xx!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx, ny, nz = size(tensor)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices(tensor)[index])
        sum = 0.0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
                x = (i-1 + p*nx)*dx
                y = (j-1 + q*ny)*dy
                z = (k-1 + s*nz)*dz
                sum += demag_tensor_xx_gpu(x,y,z,dx,dy,dz)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_tensors_kernel_yy!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx, ny, nz = size(tensor)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices(tensor)[index])
        sum = 0.0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
                x = (i-1 + p*nx)*dx
                y = (j-1 + q*ny)*dy
                z = (k-1 + s*nz)*dz
                sum += demag_tensor_yy_gpu(x,y,z,dx,dy,dz)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_tensors_kernel_zz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx, ny, nz = size(tensor)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices(tensor)[index])
        sum = 0.0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
                x = (i-1 + p*nx)*dx
                y = (j-1 + q*ny)*dy
                z = (k-1 + s*nz)*dz
                sum += demag_tensor_zz_gpu(x,y,z,dx,dy,dz)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_tensors_kernel_xy!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx, ny, nz = size(tensor)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices(tensor)[index])
        sum = 0.0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
                x = (i-1 + p*nx)*dx
                y = (j-1 + q*ny)*dy
                z = (k-1 + s*nz)*dz
                sum += demag_tensor_xy_gpu(x,y,z,dx,dy,dz)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_tensors_kernel_xz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx, ny, nz = size(tensor)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices(tensor)[index])
        sum = 0.0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
                x = (i-1 + p*nx)*dx
                y = (j-1 + q*ny)*dy
                z = (k-1 + s*nz)*dz
                sum += demag_tensor_xz_gpu(x,y,z,dx,dy,dz)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_tensors_kernel_yz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                 Nx::Int64, Ny::Int64, Nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx, ny, nz = size(tensor)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices(tensor)[index])
        sum = 0.0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
                x = (i-1 + p*nx)*dx
                y = (j-1 + q*ny)*dy
                z = (k-1 + s*nz)*dz
                sum += demag_tensor_yz_gpu(x,y,z,dx,dy,dz)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function fill_tensors_kernel!(long_tensor, tensor, tx::Bool, ty::Bool, tz::Bool)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    lnx, lny, lnz = size(long_tensor)
    nx, ny, nz = size(tensor)
    if 0 < index <= lnx*lny*lnz
        i,j,k = Tuple(CUDA.CartesianIndices(long_tensor)[index])
        if (lnx%2 == 0 && i == nx+1) || (lny%2 == 0 && j == ny+1) || (lnz%2 == 0 && k == nz+1)
          return nothing
        end
        x = (i<=nx) ? i : lnx - i + 2
        y = (j<=ny) ? j : lny - j + 2
        z = (k<=nz) ? k : lnz - k + 2
        sx = tx && (i>nx) ? -1 : 1
        sy = ty && (j>ny) ? -1 : 1
        sz = tz && (k>nz) ? -1 : 1
        @inbounds long_tensor[i,j,k] = sx*sy*sz*tensor[x,y,z]
    end
    return nothing
end

function distribute_m_kernel!(m_gpu, mx_gpu, my_gpu, mz_gpu, Ms, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(mx_gpu)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
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
        i,j,k = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        mx_gpu[i,j,k] = m_gpu[3*index-2]
        my_gpu[i,j,k] = m_gpu[3*index-1]
        mz_gpu[i,j,k] = m_gpu[3*index]
    end
    return nothing
end

function collect_h_kernel!(h, hx, hy, hz, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0< index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
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
