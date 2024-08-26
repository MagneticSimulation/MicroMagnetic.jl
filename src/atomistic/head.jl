abstract type Interaction end
"""
    AtomisticSim{T<:AbstractFloat} <: AbstractSim
"""
mutable struct AtomisticSim{T<:AbstractFloat} <: AbstractSim
    time::Float64
    mesh::AtomisticMesh
    driver::Driver
    saver::DataSaver
    spin::AbstractArray{T,1}
    prespin::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    mu_s::AbstractArray{T,1}
    pins::AbstractArray{Bool,1}
    n_total::Int64
    name::String
    driver_name::String
    interactions::Array
    save_data::Bool
    AtomisticSim{T}() where {T<:AbstractFloat} = new()
end

# HeisenbergExchange denotes the Heisenberg exchange interaction 
mutable struct HeisenbergExchange{T<:AbstractFloat} <: Interaction
    Js1::AbstractArray{T,1}  #The length of Js should be equal to the number of neighbours
    Js2::AbstractArray{T,1} #The exchange interaction for the next-nearest neigbours
    Js3::AbstractArray{T,1} #The exchange interaction for the next-next-nearest neigbours
    Js4::AbstractArray{T,1} #The exchange interaction for the next-next-next-nearest neigbours
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

# HeisenbergDMI denotes the DM interaction that can be used in 
# cubic and triangular meshes. 
mutable struct HeisenbergDMI{T<:AbstractFloat} <: Interaction
    Dij::AbstractArray{T,2}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

# DMI in Canted AFM
# HeisenbergCantedDMI denotes the DM interaction that can be used in 
# cubic and triangular meshes. 
mutable struct HeisenbergCantedDMI{T<:AbstractFloat} <: Interaction
    Dij::AbstractArray{T,2}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

# HeisenbergTubeBulkDMI denotes the Bulk DM interaction that can be used in 
# cylindrical tube mesh
mutable struct HeisenbergTubeBulkDMI{T<:AbstractFloat} <: Interaction
    D::T
    Dij::AbstractArray{T,2} # store the DM vector of a ring, which is D*r_ij
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct KagomeAnisotropy{T<:AbstractFloat} <: Interaction
    Ku::T
    ax1::Tuple
    ax2::Tuple
    ax3::Tuple
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct TubeAnisotropy{T<:AbstractFloat}
    Ku::T
    axes::AbstractArray{T,2}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct MagnetoelectricLaser{T<:AbstractFloat} <: Interaction
    lambda::T # Magnetoelectric coupling strength
    E::T #The ampitude of the electric field
    B::T #The ampitude of the magnetic field
    omega::T #the frequency of the laser
    delta::T #laser polarization is determined by delta. RCP(δ = 0),  linearly polarized (δ=π/2), LCP (δ=π)
    direction::Int64  #The direction of static field, 001, 110 or 111 
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end
