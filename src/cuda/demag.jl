using CuArrays.CUFFT
using FFTW
using LinearAlgebra

@info "Running CUFFT $(CUFFT.version())"

mutable struct DemagGPU{T<:AbstractFloat}
  nx_fft::Int64
  ny_fft::Int64
  nz_fft::Int64
  tensor_xx::CuArray{T, 3}
  tensor_yy::CuArray{T, 3}
  tensor_zz::CuArray{T, 3}
  tensor_xy::CuArray{T, 3}
  tensor_xz::CuArray{T, 3}
  tensor_yz::CuArray{T, 3}
  mx::CuArray{T, 3} #input for FFT
  my::CuArray{T, 3}
  mz::CuArray{T, 3}
  Mx::CuArray{Complex{T}, 3} #output for FFT
  My::CuArray{Complex{T}, 3}
  Mz::CuArray{Complex{T}, 3}
  Hx::CuArray{Complex{T}, 3}
  Hy::CuArray{Complex{T}, 3}
  Hz::CuArray{Complex{T}, 3}
  m_plan::Any
  h_plan::Any
  field::Array{T, 1}
  energy::Array{T, 1}
  name::String
end

function init_demag_gpu(sim::MicroSimGPU)
  mesh = sim.mesh
  max_size = max(mesh.dx, mesh.dy, mesh.dz)
  dx = Float64(mesh.dx/max_size)
  dy = Float64(mesh.dy/max_size)
  dz = Float64(mesh.dz/max_size)

  nx = mesh.nx
  ny = mesh.ny
  nz = mesh.nz

  cn = 10
  nx_fft = mesh.nx > cn ? 2*mesh.nx : 2*mesh.nx - 1
  ny_fft = mesh.ny > cn ? 2*mesh.ny : 2*mesh.ny - 1
  nz_fft = mesh.nz > cn ? 2*mesh.nz : 2*mesh.nz - 1

  Float = _cuda_using_double.x ? Float64 : Float32

  tensor_xx = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  tensor_yy = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  tensor_zz = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  tensor_xy = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  tensor_xz = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  tensor_yz = cuzeros(Float, nx_fft, ny_fft, nz_fft)

  blk, thr = CuArrays.cudims(tensor_xx)
  @cuda blocks=blk threads=thr compute_tensors_kernel!(tensor_xx, tensor_yy, tensor_zz,
                               tensor_xy, tensor_xz, tensor_yz, nx, ny, nz, dx, dy, dz)
  synchronize()
  plan = plan_rfft(tensor_xx)
  tensor_xx = real(plan*tensor_xx)
  tensor_yy = real(plan*tensor_yy)
  tensor_zz = real(plan*tensor_zz)
  tensor_xy = real(plan*tensor_xy)
  tensor_xz = real(plan*tensor_xz)
  tensor_yz = real(plan*tensor_yz)

  mx_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  my_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  mz_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  lenx = (nx_fft%2>0) ? nx : nx+1
  Mx = cuzeros(Complex{Float}, lenx, ny_fft, nz_fft)
  My = cuzeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Mz = cuzeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Hx = cuzeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Hy = cuzeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Hz = cuzeros(Complex{Float}, lenx, ny_fft, nz_fft)

  m_plan = plan_rfft(mx_gpu)
  h_plan = plan_irfft(Hx, nx_fft)

  field = zeros(Float, 3*sim.nxyz)
  energy = zeros(Float, sim.nxyz)
  demag = DemagGPU(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz,
                tensor_xy, tensor_xz, tensor_yz, mx_gpu, my_gpu, mz_gpu,
                Mx, My, Mz, Hx, Hy, Hz,
                m_plan, h_plan, field, energy, "DemagGPU")
  return demag
end

function compute_tensors_kernel!(tensor_xx, tensor_yy, tensor_zz,
                                 tensor_xy, tensor_xz, tensor_yz,
                                 nx::Int64, ny::Int64, nz::Int64,
                                 dx::Float64, dy::Float64, dz::Float64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(tensor_xx)
    if 0 < index <= nx_fft*ny_fft*nz_fft
        i,j,k = Tuple(CuArrays.CartesianIndices(tensor_xx)[index])
        if (nx_fft%2 == 0 && i == nx+1) || (ny_fft%2 == 0 && j == ny+1) || (nz_fft%2 == 0 && k == nz+1)
          return nothing
        end
        x = (i<=nx) ? (i-1)*dx : (i-nx_fft-1)*dx
        y = (j<=ny) ? (j-1)*dy : (j-ny_fft-1)*dy
        z = (k<=nz) ? (k-1)*dz : (k-nz_fft-1)*dz
        @inbounds tensor_xx[i,j,k] = demag_tensor_xx_gpu(x,y,z,dx,dy,dz)
        @inbounds tensor_yy[i,j,k] = demag_tensor_xx_gpu(y,x,z,dy,dx,dz)
        @inbounds tensor_zz[i,j,k] = demag_tensor_xx_gpu(z,y,x,dz,dy,dx)
        @inbounds tensor_xy[i,j,k] = demag_tensor_xy_gpu(x,y,z,dx,dy,dz)
        @inbounds tensor_xz[i,j,k] = demag_tensor_xy_gpu(x,z,y,dx,dz,dy)
        @inbounds tensor_yz[i,j,k] = demag_tensor_xy_gpu(y,z,x,dy,dz,dx)
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

function collect_h_kernel!(h, hx, hy, hz, nx::Int64,ny::Int64,nz::Int64)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
	mu0 = 4*pi*1e-7
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

function effective_field(demag::DemagGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

  blocks_n, threads_n = CuArrays.cudims(demag.mx)
  @cuda blocks = blocks_n threads=threads_n distribute_m_kernel!(spin, demag.mx, demag.my, demag.mz, sim.Ms, nx, ny, nz)

  #synchronize()
  mul!(demag.Mx, demag.m_plan, demag.mx)
  mul!(demag.My, demag.m_plan, demag.my)
  mul!(demag.Mz, demag.m_plan, demag.mz)

  demag.Hx .= demag.tensor_xx.*demag.Mx .+ demag.tensor_xy.*demag.My .+  demag.tensor_xz.*demag.Mz
  demag.Hy .= demag.tensor_xy.*demag.Mx .+ demag.tensor_yy.*demag.My .+  demag.tensor_yz.*demag.Mz
  demag.Hz .= demag.tensor_xz.*demag.Mx .+ demag.tensor_yz.*demag.My .+  demag.tensor_zz.*demag.Mz
  #synchronize()

  mul!(demag.mx, demag.h_plan, demag.Hx)
  mul!(demag.my, demag.h_plan, demag.Hy)
  mul!(demag.mz, demag.h_plan, demag.Hz)

  @cuda blocks = blocks_n threads=threads_n collect_h_kernel!(sim.field, demag.mx, demag.my, demag.mz, nx, ny, nz)

  blocks_n, threads_n = CuArrays.cudims(sim.nxyz)
  @cuda blocks=blocks_n threads=threads_n compute_energy_kernel!(sim.energy, sim.field, sim.field, sim.Ms, mesh.volume, sim.nxyz)
end


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
