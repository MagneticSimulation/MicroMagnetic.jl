abstract type AbstractSimGPU <:AbstractSim end
abstract type MeshGPU <: Mesh end
abstract type DriverGPU end
abstract type EnergyGPU end
abstract type MicroEnergyGPU end

mutable struct MicroSimGPU{T<:AbstractFloat} <:AbstractSimGPU
  mesh::MeshGPU
  driver::DriverGPU
  saver::DataSaver
  spin::CuArray{T, 1}
  prespin::CuArray{T, 1}
  field::CuArray{T, 1}
  energy::CuArray{T, 1}
  Ms::CuArray{T, 1}
  pins::CuArray{Bool, 1}
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

mutable struct ExchangeAnistropyGPU{T<:AbstractFloat} <: MicroEnergyGPU
   kea::CuArray{T, 1}
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

mutable struct SpatialInterfacialDMIGPU{T<:AbstractFloat} <: MicroEnergyGPU
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

mutable struct StochasticFieldGPU{F<:AbstractFloat} <: MicroEnergyGPU
    T::CuArray{F, 1}
    eta::CuArray{F, 1}
    field::Array{F, 1}
    energy::Array{F, 1}
    total_energy::F
    nsteps::Int64
    name::String
end

mutable struct TimeZeemanGPU{T<:AbstractFloat} <: MicroEnergyGPU
   time_fun::Function
   init_field::CuArray{T, 1}
   cufield::CuArray{T, 1}
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
   axis::Array{T, 1}
   Kc::T
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end
