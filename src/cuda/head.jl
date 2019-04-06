using CuArrays

abstract type MeshGPU <: Mesh end
abstract type AbstractSimGPU <:AbstractSim end
abstract type DriverGPU end
abstract type MicroEnergyGPU end

struct FDMeshGPU <: MeshGPU
  dx::FloatGPU
  dy::FloatGPU
  dz::FloatGPU
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
  volume::FloatGPU
end

function FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc=(false, false, false))
  nxyz = nx*ny*nz
  volume = dx*dy*dz
  return FDMeshGPU(FloatGPU(dx), FloatGPU(dy), FloatGPU(dz), nx, ny, nz, nxyz, pbc[1], pbc[2], pbc[3], FloatGPU(volume))
end

mutable struct Dopri5GPU
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::CuArray{FloatGPU, 1}
   omega_t::CuArray{FloatGPU, 1}
   dw_dt::CuArray{FloatGPU, 1}
   k1::CuArray{FloatGPU, 1}
   k2::CuArray{FloatGPU, 1}
   k3::CuArray{FloatGPU, 1}
   k4::CuArray{FloatGPU, 1}
   k5::CuArray{FloatGPU, 1}
   k6::CuArray{FloatGPU, 1}
   k7::CuArray{FloatGPU, 1}
   rhs_fun::Function
   succeed::Bool
end

mutable struct LLG_GPU <: DriverGPU
  precession::Bool
  alpha::Float64
  gamma::Float64
  ode::Dopri5GPU
  tol::Float64
  field::CuArray{FloatGPU, 1}
end

mutable struct MicroSimGPU <:AbstractSimGPU
  mesh::FDMeshGPU
  driver::DriverGPU
  saver::DataSaver
  spin::CuArray{FloatGPU, 1}
  prespin::CuArray{FloatGPU, 1}
  field::CuArray{FloatGPU, 1}
  energy::CuArray{FloatGPU, 1}
  Ms::CuArray{FloatGPU, 1}
  nxyz::Int64
  blocks::Int64
  threads::Int64
  name::String
  interactions::Array{Any, 1}
end

mutable struct ExchangeGPU <: MicroEnergyGPU
   A::Float64
   field::Array{FloatGPU, 1}
   energy::Array{FloatGPU, 1}
   total_energy::FloatGPU
   name::String
end

mutable struct ZeemanGPU <: MicroEnergyGPU
   Hx::Float64
   Hy::Float64
   Hz::Float64
   field::Array{FloatGPU, 1}
   field_gpu::CuArray{FloatGPU, 1}
   energy::Array{FloatGPU, 1}
   total_energy::FloatGPU
   name::String
end

mutable struct AnisotropyGPU <: MicroEnergyGPU
   Ku::CuArray{FloatGPU, 1}
   axis::Tuple
   field::Array{FloatGPU, 1}
   energy::Array{FloatGPU, 1}
   total_energy::FloatGPU
   name::String
end
