abstract type InteractionGPU end

#AtomicSimGPU is the overall struct for atomistic simulation.
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
  n_total::Int64
  blocks::Int64
  threads::Int64
  name::String
  driver_name::String
  interactions::Array{Any, 1}
  save_data::Bool
  AtomicSimGPU{T}() where {T<:AbstractFloat} = new()
end

# HeisenbergExchange denotes the Heisenberg exchange interaction 
mutable struct HeisenbergExchange{T<:AbstractFloat} <: InteractionGPU
   Js::CuArray{T, 1}  #The length of Js should be equal to the number of neighbours
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end


mutable struct NextHeisenbergExchange{T<:AbstractFloat} <: InteractionGPU
   Js::CuArray{T, 1}  #The length of Js should be equal to the number of neighbours
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct NextNextHeisenbergExchange{T<:AbstractFloat} <: InteractionGPU
   Js::CuArray{T, 1}  #The length of Js should be equal to the number of neighbours
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

mutable struct NextNextNextHeisenbergExchange{T<:AbstractFloat} <: InteractionGPU
   Js::CuArray{T,1}  #The length of Js should be equal to the number of neighbours
   field::Array{T,1}
   energy::Array{T,1}
   total_energy::T
   name::String
end
# HeisenbergBulkDMI denotes the Bulk DM interaction that can be used in 
# cubic and triangular meshes. 
mutable struct HeisenbergBulkDMI{T<:AbstractFloat} <: InteractionGPU
   D::T 
   field::Array{T, 1}
   energy::Array{T, 1}
   total_energy::T
   name::String
end

# HeisenbergTubeBulkDMI denotes the Bulk DM interaction that can be used in 
# cylindrical tube mesh
mutable struct HeisenbergTubeBulkDMI{T<:AbstractFloat} <: InteractionGPU
   D::T
   Dij::CuArray{T, 2} # store the DM vector of a ring, which is D*r_ij
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

mutable struct TubeAnisotropy{T<:AbstractFloat}
   Ku::T
   axes::CuArray{T, 2}
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