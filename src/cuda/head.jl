abstract type MeshGPU <: Mesh end
abstract type AbstractSimGPU <:AbstractSim end
abstract type DriverGPU end
abstract type MicroEnergyGPU end

struct FDMeshGPU{T <: AbstractFloat} <: MeshGPU
  dx::T
  dy::T
  dz::T
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
  volume::T
end

mutable struct Dopri5GPU{T<:AbstractFloat} <: AbstractDopri5
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::CuArray{T, 1}
   omega_t::CuArray{T, 1}
   dw_dt::CuArray{T, 1}
   k1::CuArray{T, 1}
   k2::CuArray{T, 1}
   k3::CuArray{T, 1}
   k4::CuArray{T, 1}
   k5::CuArray{T, 1}
   k6::CuArray{T, 1}
   k7::CuArray{T, 1}
   rhs_fun::Function
   succeed::Bool
end

mutable struct MicroSimGPU{T<:AbstractFloat} <:AbstractSimGPU
  mesh::FDMeshGPU
  driver::DriverGPU
  saver::DataSaver
  spin::CuArray{T, 1}
  prespin::CuArray{T, 1}
  field::CuArray{T, 1}
  energy::CuArray{T, 1}
  Ms::CuArray{T, 1}
  total_energy::T
  nxyz::Int64
  blocks::Int64
  threads::Int64
  name::String
  interactions::Array{Any, 1}
end

mutable struct ExchangeGPU{T<:AbstractFloat} <: MicroEnergyGPU
   A::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct BulkDMIGPU{T<:AbstractFloat} <: MicroEnergyGPU
   D::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct ZeemanGPU{T<:AbstractFloat} <: MicroEnergyGPU
   Hx::Float64
   Hy::Float64
   Hz::Float64
   field::Array{T, 1}
   field_gpu::CuArray{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct AnisotropyGPU{T<:AbstractFloat} <: MicroEnergyGPU
   Ku::CuArray{T, 1}
   axis::Tuple
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end
