abstract type InteractionGPU end

mutable struct AtomicSimGPU{T<:AbstractFloat} <:AbstractSimGPU
  mesh::MeshGPU
  driver::DriverGPU
  saver::DataSaver
  spin::CuArray{T, 1}
  prespin::CuArray{T, 1}
  field::CuArray{T, 1}
  energy::CuArray{T, 1}
  mu_s::CuArray{T, 1}
  pins::CuArray{Bool, 1}
  total_energy::T
  nxyz::Int64
  blocks::Int64
  threads::Int64
  name::String
  driver_name::String
  interactions::Array{Any, 1}
  save_data::Bool
  AtomicSimGPU{T}() where {T<:AbstractFloat} = new()
end

mutable struct HeisenbergExchange{T<:AbstractFloat} <: InteractionGPU
   Js::CuArray{T, 1}  #The length of Js should be equal to the number of neighbours
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct HeisenbergBulkDMI{T<:AbstractFloat} <: InteractionGPU
   D::T 
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct KagomeAnisotropy{T<:AbstractFloat} <: InteractionGPU
   Ku::T
   ax1::Tuple
   ax2::Tuple
   ax3::Tuple
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end


mutable struct MagnetoelectricLaser{T<:AbstractFloat} <: InteractionGPU
   lambda::T # Magnetoelectric coupling strength
   E::T #The ampitude of the electric field
   B::T #The ampitude of the magnetic field
   omega:: T #the frequency of the laser
   delta:: T #laser polarization is determined by delta. RCP(δ = 0),  linearly polarized (δ=π/2), LCP (δ=π)
   direction::Int64  #The direction of static field, 001, 110 or 111 
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end