abstract type Driver end
abstract type AbstractSim end
abstract type MicroEnergy end
abstract type Integrator end
abstract type IntegratorCayley <: Integrator end

"""
    NumberOrArrayOrFunction

In Micromagnetics, typical parameters such as saturation magnetization `Ms` and exchange stiffness constant `A` are constant.
However, there are many cases that a spatial `Ms` is needed. For instance, if the simulated system is a circular disk the Ms
in the corners should be set to zero. If the simulated system contains mutiple materials, the exchange constant `A` should be
spatial as well. The Union `NumberOrArrayOrFunction` is designed to deal with such situations. As indicated from its name, it
means that the input parameter could be a number or an array or a function:
  - Number: should be Real.
  - Array: the length of the array should be `N` where `N` is the total spin number of the system.
  - Function: the parameter of the function should be `(i,j,k,dx,dy,dz)` where `i,j,k` is the cell index and `dx,dy,dz` is the cellsize.
    The return value of the function should be a real number. For example,

    ```julia
    function circular_Ms(i,j,k,dx,dy,dz)
        if (i-50.5)^2 + (j-50.5)^2 <= 50^2
            return 8.6e5
        end
        return 0.0
    end
    ```
"""
NumberOrArrayOrFunction = Union{Number,Array,Function}

"""
    NumberOrTupleOrArrayOrFunction

In micromagnetics, there are also cases where the input parameters can be either scalars or vectors and vary with space. For example,
the parameters for the DMI could be a const for bulk DMI or interfacial DMI. In some materials, the DMI const may differ in different 
directions and thus a tuple with three numbers is required. In MicroMagnetic, the union `NumberOrTupleOrArrayOrFunction` is designed to deal 
with such situations. Similar to `NumberOrArrayOrFunction`, `NumberOrTupleOrArrayOrFunction` means that the input parameter could be 
a number, a tuple, an array or a function:
  - Number: should be Real.
  - Tuple: should be Real with length 3. For example, `(1,2e-5,0)`.
  - Array: the length of the array should be `N` or `3N` where `N` is the total spin number of the system.
  - Function: the parameter of the function should be `(i,j,k,dx,dy,dz)` and the return value should be a tuple with length 3.

    For example,
    ```julia
    function uniform_m0(i,j,k,dx,dy,dz)
        return (0,0,1)
    end
    ```
"""
NumberOrTupleOrArrayOrFunction = Union{Number,Tuple,Array,Function}

"""
    NumberOrArray

Similar to Union `NumberOrArrayOrFunction`, the Union `NumberOrArray` is designed to deal with cases that a number
or an array is needed.
"""
NumberOrArray = Union{Number,Array}

"""
    ArrayOrFunction

Similar to Union `NumberOrArrayOrFunction`, the Union `ArrayOrFunction` is designed to deal with cases that
an array or a function is needed.
"""
ArrayOrFunction = Union{Array,Function}

"""
    TupleOrArrayOrFunction

Similar to `NumberOrArrayOrFunction`, `TupleOrArrayOrFunction` means that the input parameter could be a tuple or
an array or a function:
  - Tuple: should be Real with length 3. For example, `(0,0,1e5)`.
  - Array: the length of the array should be `3N` where `N` is the total spin number of the system.
  - Function: the parameter of the function should be `(i,j,k,dx,dy,dz)` and the return value should be a tuple with length 3.
    For example,
    ```julia
    function uniform_m0(i,j,k,dx,dy,dz)
        return (0,0,1)
    end
    ```
"""
TupleOrArrayOrFunction = Union{Tuple,Array,Function}

"""
    TupleOrArray

The Union `TupleOrArray` is designed to deal with cases that a Tuple or an array is needed.
"""
TupleOrArray = Union{Tuple,Array}

mutable struct DataSaver
    name::String
    header_saved::Bool
    t::Float64
    nsteps::Int
    items::Array
end

mutable struct SaverItem
    name::Union{String,Tuple}
    unit::Union{String,Tuple}
    result::Function
end

"""
    MicroSim{T<:AbstractFloat} <: AbstractSim
"""
mutable struct MicroSim{T<:AbstractFloat} <: AbstractSim
    time::Float64
    mesh::FDMesh
    driver::Driver
    saver::DataSaver
    spin::AbstractArray{T,1}
    prespin::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    mu0_Ms::AbstractArray{T,1}
    pins::AbstractArray{Bool,1}
    n_total::Int64
    name::String
    driver_name::String
    interactions::Array{MicroEnergy}
    save_data::Bool
    MicroSim{T}() where {T<:AbstractFloat} = new()
end

mutable struct Zeeman{T<:AbstractFloat} <: MicroEnergy
    H0::Tuple
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct TimeZeeman{T<:AbstractFloat} <: MicroEnergy
    time_fun::Function
    init_field::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct Anisotropy{T<:AbstractFloat} <: MicroEnergy
    Ku::AbstractArray{T,1}
    axis::Tuple{T,T,T}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct CubicAnisotropy{T<:AbstractFloat} <: MicroEnergy
    Kc::AbstractArray{T,1}
    axis1::Tuple
    axis2::Tuple
    axis3::Tuple
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct SpatialExchange{T<:AbstractFloat} <: MicroEnergy
    A::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct UniformExchange{T<:AbstractFloat} <: MicroEnergy
    Ax::Real
    Ay::Real
    Az::Real
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct InterlayerExchange{T<:AbstractFloat} <: MicroEnergy
    J::Float64
    k1::Int32
    k2::Int32
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct InterlayerDMI{T<:AbstractFloat} <: MicroEnergy
    Dx::T
    Dy::T
    Dz::T
    k1::Int32
    k2::Int32
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct BulkDMI{T<:AbstractFloat} <: MicroEnergy
    Dx::Real
    Dy::Real
    Dz::Real
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct SpatialBulkDMI{T<:AbstractFloat} <: MicroEnergy
    D::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct InterfacialDMI{T<:AbstractFloat} <: MicroEnergy
    D::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct StochasticField{T<:AbstractFloat} <: MicroEnergy
    temperature::AbstractArray{T,1}
    eta::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    nsteps::Int64
    name::String
    k_B::Float64
    scaling_fun::Function
    average_temperature::T
    scaling_factor::T
end
