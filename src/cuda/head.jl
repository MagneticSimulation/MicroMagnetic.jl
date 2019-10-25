abstract type AbstractSimGPU <:AbstractSim end
abstract type DriverGPU end
abstract type MicroEnergyGPU end

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
  save_data::Bool
  MicroSimGPU{T}() where {T<:AbstractFloat} = new()
end


mutable struct ExchangeGPU{T<:AbstractFloat} <: MicroEnergyGPU
   A::CuArray{T, 1}
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct Vector_ExchangeGPU{T<:AbstractFloat} <: MicroEnergyGPU
   A::CuArray{T, 1}
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct ExchangeRKKYGPU{T<:AbstractFloat} <: MicroEnergyGPU
   sigma::T
   Delta::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct BulkDMIGPU{T<:AbstractFloat} <: MicroEnergyGPU
   Dx::T
   Dy::T
   Dz::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct InterfacialDMIGPU{T<:AbstractFloat} <: MicroEnergyGPU
   D::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct SpatialBulkDMIGPU{T<:AbstractFloat} <: MicroEnergyGPU
   D::CuArray{T, 1}
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct ZeemanGPU{T<:AbstractFloat} <: MicroEnergyGPU
   field::Array{T, 1}
   energy::Array{T, 1}
   cufield::CuArray{T, 1}
   total_energy::T
   name::String
end

mutable struct TimeZeemanGPU{T<:AbstractFloat} <: MicroEnergyGPU
   time_fun::Function
   init_field::CuArray{T, 1}
   field::Array{T, 1}
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

mutable struct CubicAnisotropyGPU{T<:AbstractFloat} <: MicroEnergyGPU
   Kc::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end
