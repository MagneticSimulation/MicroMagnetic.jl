using CuArrays, CUDAnative, CUDAdrv
#CuArrays.allowscalar(false)
using CuArrays.CUFFT
using FFTW
using LinearAlgebra

@info "Running CUFFT $(CUFFT.version())"

mutable struct DemagGPU
  nx_fft::Int64
  ny_fft::Int64
  nz_fft::Int64
  tensor_xx::CuArray{Float64, 3}
  tensor_yy::CuArray{Float64, 3}
  tensor_zz::CuArray{Float64, 3}
  tensor_xy::CuArray{Float64, 3}
  tensor_xz::CuArray{Float64, 3}
  tensor_yz::CuArray{Float64, 3}
  m_cpu::Array{Float64, 1}
  m_gpu::CuArray{Float64, 1}  #receive spin from CPU
  mx_gpu::CuArray{Float64, 3} #input for FFT
  my_gpu::CuArray{Float64, 3}
  mz_gpu::CuArray{Float64, 3}
  Mx::CuArray{Complex{Float64}, 3} #output for FFT
  My::CuArray{Complex{Float64}, 3}
  Mz::CuArray{Complex{Float64}, 3}
  Hx::CuArray{Complex{Float64}, 3}
  Hy::CuArray{Complex{Float64}, 3}
  Hz::CuArray{Complex{Float64}, 3}
  m_plan::Any
  h_plan::Any
  field::Array{Float64, 1}
  energy::Array{Float64, 1}
  name::String
end

function init_demag_gpu(sim::SimData)
  mesh = sim.mesh
  max_size = max(mesh.dx, mesh.dy, mesh.dz)
  dx = mesh.dx/max_size
  dy = mesh.dy/max_size
  dz = mesh.dz/max_size

  nx = mesh.nx
  ny = mesh.ny
  nz = mesh.nz

  cn = 10
  nx_fft = mesh.nx > cn ? 2*mesh.nx : 2*mesh.nx - 1
  ny_fft = mesh.ny > cn ? 2*mesh.ny : 2*mesh.ny - 1
  nz_fft = mesh.nz > cn ? 2*mesh.nz : 2*mesh.nz - 1

  tensor_xx = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  tensor_yy = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  tensor_zz = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  tensor_xy = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  tensor_xz = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  tensor_yz = cuzeros(Float64, nx_fft, ny_fft, nz_fft)

  numblocks, numthreads = CuArrays.cudims(tensor_xx)
  @cuda blocks = numblocks threads=numthreads compute_tensors_kernel!(tensor_xx, tensor_yy, tensor_zz,
                                            tensor_xy, tensor_xz, tensor_yz, nx, ny, nz, dx, dy, dz)
  synchronize()
  plan = plan_rfft(tensor_xx)
  tensor_xx = real(plan*tensor_xx)
  tensor_yy = real(plan*tensor_yy)
  tensor_zz = real(plan*tensor_zz)
  tensor_xy = real(plan*tensor_xy)
  tensor_xz = real(plan*tensor_xz)
  tensor_yz = real(plan*tensor_yz)

  m_cpu = zeros(Float64, 3*sim.nxyz)
  m_gpu = cuzeros(Float64, 3*sim.nxyz)
  mx_gpu = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  my_gpu = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  mz_gpu = cuzeros(Float64, nx_fft, ny_fft, nz_fft)
  lenx = (nx_fft%2>0) ? nx : nx+1
  Mx = cuzeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  My = cuzeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Mz = cuzeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Hx = cuzeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Hy = cuzeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Hz = cuzeros(Complex{Float64}, lenx, ny_fft, nz_fft)

  m_plan = plan_rfft(mx_gpu)
  h_plan = plan_irfft(Hx, nx_fft)

  field = zeros(Float64, 3*sim.nxyz)
  energy = zeros(Float64, sim.nxyz)
  demag = DemagGPU(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz,
                tensor_xy, tensor_xz, tensor_yz, m_cpu, m_gpu, mx_gpu, my_gpu, mz_gpu,
                Mx, My, Mz, Hx, Hy, Hz,
                m_plan, h_plan, field, energy, "DemagGPU")
  return demag
end


function compute_tensors_kernel!(tensor_xx, tensor_yy, tensor_zz, tensor_xy, tensor_xz, tensor_yz, nx, ny, nz, dx, dy, dz)
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
        tensor_xx[i,j,k] = demag_tensor_xx_gpu(x,y,z,dx,dy,dz)
        tensor_yy[i,j,k] = demag_tensor_xx_gpu(y,x,z,dy,dx,dz)
        tensor_zz[i,j,k] = demag_tensor_xx_gpu(z,y,x,dz,dy,dx)
        tensor_xy[i,j,k] = demag_tensor_xy_gpu(x,y,z,dx,dy,dz)
        tensor_xz[i,j,k] = demag_tensor_xy_gpu(x,z,y,dx,dz,dy)
        tensor_yz[i,j,k] = demag_tensor_xy_gpu(y,z,x,dy,dz,dx)
    end
    return nothing
end

function distribute_m_kernel!(m_gpu, mx_gpu, my_gpu, mz_gpu, nx,ny,nz)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(mx_gpu)
    if 0< index <= nx*ny*nz
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
        mx_gpu[i,j,k] = m_gpu[3*index-2]
        my_gpu[i,j,k] = m_gpu[3*index-1]
        mz_gpu[i,j,k] = m_gpu[3*index]
    end
    return nothing
end


function collect_h_kernel!(m_gpu, mx_gpu, my_gpu, mz_gpu, nx,ny,nz)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    nx_fft, ny_fft, nz_fft = size(mx_gpu)
    if 0< index <= nx*ny*nz
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
        m_gpu[3*index-2] = -1.0*mx_gpu[i,j,k]
        m_gpu[3*index-1] = -1.0*my_gpu[i,j,k]
        m_gpu[3*index] = -1.0*mz_gpu[i,j,k]
    end
    return nothing
end

function effective_field(demag::DemagGPU, sim::SimData, spin::Array{Float64, 1}, t::Float64)
  nx, ny, nz = sim.mesh.nx, sim.mesh.ny, sim.mesh.nz
  for i = 1:sim.nxyz
      j = 3*i
      demag.m_cpu[j-2]  = spin[j-2]*sim.Ms[i]
      demag.m_cpu[j-1]  = spin[j-1]*sim.Ms[i]
      demag.m_cpu[j]  = spin[j]*sim.Ms[i]
  end
  copyto!(demag.m_gpu, demag.m_cpu);
  blocks_n, threads_n = CuArrays.cudims(demag.mx_gpu)
  @cuda blocks = blocks_n threads=threads_n distribute_m_kernel!(demag.m_gpu, demag.mx_gpu, demag.my_gpu, demag.mz_gpu, nx, ny, nz)
  synchronize()
  mul!(demag.Mx, demag.m_plan, demag.mx_gpu)
  mul!(demag.My, demag.m_plan, demag.my_gpu)
  mul!(demag.Mz, demag.m_plan, demag.mz_gpu)


  demag.Hx .= demag.tensor_xx.*demag.Mx .+ demag.tensor_xy.*demag.My .+  demag.tensor_xz.*demag.Mz
  demag.Hy .= demag.tensor_xy.*demag.Mx .+ demag.tensor_yy.*demag.My .+  demag.tensor_yz.*demag.Mz
  demag.Hz .= demag.tensor_xz.*demag.Mx .+ demag.tensor_yz.*demag.My .+  demag.tensor_zz.*demag.Mz
  synchronize()

  mul!(demag.mx_gpu, demag.h_plan, demag.Hx)
  mul!(demag.my_gpu, demag.h_plan, demag.Hy)
  mul!(demag.mz_gpu, demag.h_plan, demag.Hz)

  @cuda blocks = blocks_n threads=threads_n collect_h_kernel!(demag.m_gpu, demag.mx_gpu, demag.my_gpu, demag.mz_gpu, nx, ny, nz)
  copyto!(demag.field, demag.m_gpu);

  mu0 = 4*pi*1e-7
  volume = sim.mesh.volume
  field = demag.field
  for i=1:sim.nxyz
    j = 3*i
    demag.energy[i] = -0.5*mu0*volume*sim.Ms[i]*(field[j-2]*spin[j-2] + field[j-1]*spin[j-1] + field[j]*spin[j])
  end
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
