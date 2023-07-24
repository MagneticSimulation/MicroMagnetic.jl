using CUDA.CUFFT
using FFTW
using LinearAlgebra

function dipolar_tensor_xx_gpu(x::Float64, y::Float64, z::Float64)
    R = x*x + y*y + z*z
    if R == 0
        return 0.0
    else
        return -(2 * x*x - y*y - z*z) / (R * R * CUDA.sqrt(R))
    end
end

function dipolar_tensor_xy_gpu(x::Float64, y::Float64, z::Float64)
    R = x*x + y*y + z*z
    if R == 0
        return 0.0
    else
        return - 3*x*y / (R * R * CUDA.sqrt(R))
    end
end

function dipolar_tensor_yy_gpu(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xx_gpu(y,x,z);
end

function dipolar_tensor_zz_gpu(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xx_gpu(z,y,x);
end

function dipolar_tensor_xz_gpu(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xy_gpu(x,z,y);
end

function dipolar_tensor_yz_gpu(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xy_gpu(y,z,x)
end

function compute_dipolar_tensors_kernel_xx!(tensor, dx::Float64, dy::Float64, dz::Float64,
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
                sum += dipolar_tensor_xx_gpu(x,y,z)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_dipolar_tensors_kernel_yy!(tensor, dx::Float64, dy::Float64, dz::Float64,
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
                sum += dipolar_tensor_yy_gpu(x,y,z)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_dipolar_tensors_kernel_zz!(tensor, dx::Float64, dy::Float64, dz::Float64,
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
                sum += dipolar_tensor_zz_gpu(x,y,z)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_dipolar_tensors_kernel_xy!(tensor, dx::Float64, dy::Float64, dz::Float64,
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
                sum += dipolar_tensor_xy_gpu(x,y,z)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_dipolar_tensors_kernel_xz!(tensor, dx::Float64, dy::Float64, dz::Float64,
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
                sum += dipolar_tensor_xz_gpu(x,y,z)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function compute_dipolar_tensors_kernel_yz!(tensor, dx::Float64, dy::Float64, dz::Float64,
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
                sum += dipolar_tensor_yz_gpu(x,y,z)
        end
        @inbounds tensor[i,j,k] = sum
    end
    return nothing
end

function distribute_m_atomistic_kernel!(m_gpu, mx_gpu, my_gpu, mz_gpu, mu_s, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(mx_gpu)
    if 0 < index <= nx*ny*nz
        i,j,k = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
		p = 3*index-2
        @inbounds mx_gpu[i,j,k] = m_gpu[p]*mu_s[index]*1e20 
        @inbounds my_gpu[i,j,k] = m_gpu[p+1]*mu_s[index]*1e20
        @inbounds mz_gpu[i,j,k] = m_gpu[p+2]*mu_s[index]*1e20
    end
    return nothing
end

function compute_energy_kernel!(energy::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1},
                               mu_s::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0< index <= nxyz
        j = 3*index-2
        @inbounds energy[index] = -0.5*mu_s[index]*(m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2])
    end
    return nothing
end
