using CuArrays.CUFFT
using FFTW
using LinearAlgebra

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
  total_energy::T
  name::String
end

function init_demag_gpu(sim::MicroSimGPU, Nx::Int, Ny::Int, Nz::Int)
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

  mx_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  my_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  mz_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)

  blk, thr = CuArrays.cudims(mx_gpu)
  @cuda blocks=blk threads=thr compute_tensors_kernel_xx!(mx_gpu, my_gpu, mz_gpu,
                                                          nx, ny, nz, dx, dy, dz,
                                                          Nx, Ny, Nz)

  #synchronize()
  plan = plan_rfft(mx_gpu)
  tensor_xx = real(plan*mx_gpu)
  tensor_yy = real(plan*my_gpu)
  tensor_zz = real(plan*mz_gpu)

  @cuda blocks=blk threads=thr compute_tensors_kernel_xy!(mx_gpu, my_gpu, mz_gpu,
                                                          nx, ny, nz, dx, dy, dz,
                                                          Nx, Ny, Nz)
  #synchronize()
  tensor_xy = real(plan*mx_gpu)
  tensor_xz = real(plan*my_gpu)
  tensor_yz = real(plan*mz_gpu)

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
                m_plan, h_plan, field, energy, Float(0), "DemagGPU")
  return demag
end

mutable struct DemagGPUII{T<:AbstractFloat}
  nx_fft::Int64
  ny_fft::Int64
  nz_fft::Int64
  tensor_xx::CuArray{T, 3}
  tensor_yy::CuArray{T, 3}
  tensor_zz::CuArray{T, 3}
  tensor_xy::CuArray{T, 3}
  tensor_xz::CuArray{T, 3}
  tensor_yz::CuArray{T, 3}
  spin::CuArray{T, 1}
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

#simply make a copy for now
function init_demag_gpu_II(sim::MicroSim, Nx::Int, Ny::Int, Nz::Int)
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

  mx_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  my_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)
  mz_gpu = cuzeros(Float, nx_fft, ny_fft, nz_fft)

  blk, thr = CuArrays.cudims(mx_gpu)
  @cuda blocks=blk threads=thr compute_tensors_kernel_xx!(mx_gpu, my_gpu, mz_gpu,
                                                          nx, ny, nz, dx, dy, dz,
                                                          Nx, Ny, Nz)

  #synchronize()
  plan = plan_rfft(mx_gpu)
  tensor_xx = real(plan*mx_gpu)
  tensor_yy = real(plan*my_gpu)
  tensor_zz = real(plan*mz_gpu)

  @cuda blocks=blk threads=thr compute_tensors_kernel_xy!(mx_gpu, my_gpu, mz_gpu,
                                                          nx, ny, nz, dx, dy, dz,
                                                          Nx, Ny, Nz)
  #synchronize()
  tensor_xy = real(plan*mx_gpu)
  tensor_xz = real(plan*my_gpu)
  tensor_yz = real(plan*mz_gpu)

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
  spin = cuzeros(Float, 3*sim.nxyz)
  demag = DemagGPUII(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz,
                tensor_xy, tensor_xz, tensor_yz, spin, mx_gpu, my_gpu, mz_gpu,
                Mx, My, Mz, Hx, Hy, Hz,
                m_plan, h_plan, field, energy, "DemagGPUII")
  return demag
end

function effective_field(demag::DemagGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

  fill!(demag.mx, 0)
  fill!(demag.my, 0)
  fill!(demag.mz, 0)

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
  @cuda blocks=blocks_n threads=threads_n compute_energy_kernel!(sim.energy, spin, sim.field, sim.Ms, mesh.volume, sim.nxyz)
  return nothing
end

function effective_field(demag::DemagGPUII, sim::MicroSim, spin::Array{Float64, 1}, t::Float64) where {T<:AbstractFloat}
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz
  for i = 1:sim.nxyz
      j = 3*i-2
      demag.field[j]  = spin[j]*sim.Ms[i]
      demag.field[j+1]  = spin[j+1]*sim.Ms[i]
      demag.field[j+2]  = spin[j+2]*sim.Ms[i]
  end
  copyto!(demag.spin, demag.field)

  fill!(demag.mx, 0)
  fill!(demag.my, 0)
  fill!(demag.mz, 0)

  blocks_n, threads_n = CuArrays.cudims(demag.mx)
  @cuda blocks = blocks_n threads=threads_n distribute_m_kernel_II!(demag.spin, demag.mx, demag.my, demag.mz, nx, ny, nz)

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

  @cuda blocks = blocks_n threads=threads_n collect_h_kernel!(demag.spin, demag.mx, demag.my, demag.mz, nx, ny, nz)
  copyto!(demag.field, demag.spin);

  mu0 = 4*pi*1e-7
  volume = sim.mesh.volume
  field = demag.field
  for i=1:sim.nxyz
    j = 3*i
    demag.energy[i] = -0.5*mu0*volume*sim.Ms[i]*(field[j-2]*spin[j-2] + field[j-1]*spin[j-1] + field[j]*spin[j])
  end
  return nothing
end
