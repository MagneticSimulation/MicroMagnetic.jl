using CUDA.CUFFT
using FFTW
using LinearAlgebra

function init_demag_gpu(sim::AtomicSimGPU, Nx::Int, Ny::Int, Nz::Int)
  mesh = sim.mesh
  max_size = max(mesh.dx, mesh.dy, mesh.dz)
  #dx = Float64(mesh.dx/max_size)
  #dy = Float64(mesh.dy/max_size)
  #dz = Float64(mesh.dz/max_size)
  dx = mesh.dx*1e9
  dy = mesh.dy*1e9
  dz = mesh.dz*1e9

  nx = mesh.nx
  ny = mesh.ny
  nz = mesh.nz

  cn = 10
  nx_fft = mesh.nx > cn ? 2*mesh.nx : 2*mesh.nx - 1
  ny_fft = mesh.ny > cn ? 2*mesh.ny : 2*mesh.ny - 1
  nz_fft = mesh.nz > cn ? 2*mesh.nz : 2*mesh.nz - 1

  Float = _cuda_using_double.x ? Float64 : Float32

  mx_gpu = CUDA.zeros(Float, nx_fft, ny_fft, nz_fft)
  my_gpu = CUDA.zeros(Float, nx_fft, ny_fft, nz_fft)
  mz_gpu = CUDA.zeros(Float, nx_fft, ny_fft, nz_fft)
  plan = plan_rfft(mx_gpu)

  tensor = CUDA.zeros(Float, nx, ny, nz)

  blk1, thr1 = cudims(tensor)
  blk2, thr2 = cudims(mx_gpu)

  #Nxx
  @cuda blocks=blk1 threads=thr1 compute_dipolar_tensors_kernel_xx!(tensor, dx, dy, dz, Nx, Ny, Nz)
  @cuda blocks=blk2 threads=thr2 fill_tensors_kernel!(mx_gpu, tensor, false, false, false)
  tensor_xx = real(plan*mx_gpu)
  #Nyy
  @cuda blocks=blk1 threads=thr1 compute_dipolar_tensors_kernel_yy!(tensor, dx, dy, dz, Nx, Ny, Nz)
  @cuda blocks=blk2 threads=thr2 fill_tensors_kernel!(my_gpu, tensor, false, false, false)
  tensor_yy = real(plan*my_gpu)
  #Nzz
  @cuda blocks=blk1 threads=thr1 compute_dipolar_tensors_kernel_zz!(tensor, dx, dy, dz, Nx, Ny, Nz)
  @cuda blocks=blk2 threads=thr2 fill_tensors_kernel!(mz_gpu, tensor, false, false, false)
  tensor_zz = real(plan*mz_gpu)
  #Nxy
  @cuda blocks=blk1 threads=thr1 compute_dipolar_tensors_kernel_xy!(tensor, dx, dy, dz, Nx, Ny, Nz)
  @cuda blocks=blk2 threads=thr2 fill_tensors_kernel!(mx_gpu, tensor, true, true, false)
  tensor_xy = real(plan*mx_gpu)
  #Nxz
  @cuda blocks=blk1 threads=thr1 compute_dipolar_tensors_kernel_xz!(tensor, dx, dy, dz, Nx, Ny, Nz)
  @cuda blocks=blk2 threads=thr2 fill_tensors_kernel!(my_gpu, tensor, true, false, true)
  tensor_xz = real(plan*my_gpu)
  #Nxz
  @cuda blocks=blk1 threads=thr1 compute_dipolar_tensors_kernel_yz!(tensor, dx, dy, dz, Nx, Ny, Nz)
  @cuda blocks=blk2 threads=thr2 fill_tensors_kernel!(mz_gpu, tensor, false, true, true)
  tensor_yz = real(plan*mz_gpu)

  lenx = (nx_fft%2>0) ? nx : nx+1
  Mx = CUDA.zeros(Complex{Float}, lenx, ny_fft, nz_fft)
  My = CUDA.zeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Mz = CUDA.zeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Hx = CUDA.zeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Hy = CUDA.zeros(Complex{Float}, lenx, ny_fft, nz_fft)
  Hz = CUDA.zeros(Complex{Float}, lenx, ny_fft, nz_fft)

  m_plan = plan_rfft(mx_gpu)
  h_plan = plan_irfft(Hx, nx_fft)

  field = zeros(Float, 3*sim.n_nodes)
  energy = zeros(Float, sim.n_nodes)
  demag = DemagGPU(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz,
                tensor_xy, tensor_xz, tensor_yz, mx_gpu, my_gpu, mz_gpu,
                Mx, My, Mz, Hx, Hy, Hz,
                m_plan, h_plan, field, energy, Float(0), "DemagGPU")
  return demag
end


function effective_field(demag::DemagGPU, sim::AtomicSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

  fill!(demag.mx, 0)
  fill!(demag.my, 0)
  fill!(demag.mz, 0)

  blocks_n, threads_n = cudims(demag.mx)
  @cuda blocks = blocks_n threads=threads_n distribute_m_atomistic_kernel!(spin, demag.mx, demag.my, demag.mz, sim.mu_s, nx, ny, nz)

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

  blocks_n, threads_n = cudims(sim.n_nodes)
  @cuda blocks=blocks_n threads=threads_n compute_energy_kernel!(sim.energy, spin, sim.field, sim.mu_s, sim.n_nodes)
  return nothing
end

