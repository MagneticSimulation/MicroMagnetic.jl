abstract type Driver end
abstract type AbstractSim end
abstract type MicroEnergy end
abstract type Interaction end
abstract type Integrator end
abstract type AbstractDopri5 <:Integrator end

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
NumberOrArrayOrFunction = Union{Number, Array, Function}


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
TupleOrArrayOrFunction = Union{Tuple, Array, Function}

mutable struct Dopri5 <: AbstractDopri5
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::Array{Float64, 1}
   omega_t::Array{Float64, 1}
   dw_dt::Array{Float64, 1}
   k1::Array{Float64, 1}
   k2::Array{Float64, 1}
   k3::Array{Float64, 1}
   k4::Array{Float64, 1}
   k5::Array{Float64, 1}
   k6::Array{Float64, 1}
   k7::Array{Float64, 1}
   rhs_fun::Function
   succeed::Bool
end

mutable struct DataSaver
  name::String
  t::Float64
  nsteps::Int64
  header_saved::Bool
  headers::Array  #string or tuple<string> array
  units::Array #string or tuple<string> array
  results::Array  #function array
end

mutable struct MicroSim <: AbstractSim
  mesh::FDMesh
  driver::Driver
  saver::DataSaver
  spin::Array{Float64, 1}
  prespin::Array{Float64, 1}
  field::Array{Float64, 1}
  energy::Array{Float64, 1}
  Ms::Array{Float64, 1}
  nxyz::Int64
  name::String
  interactions::Array
  save_data::Bool
  MicroSim() = new()
end

mutable struct AtomicSim <: AbstractSim
  mesh::Mesh
  driver::Driver
  saver::DataSaver
  spin::Array{Float64, 1}
  prespin::Array{Float64, 1}
  field::Array{Float64, 1}
  energy::Array{Float64, 1}
  mu_s::Array{Float64, 1}
  nxyz::Int64
  name::String
  interactions::Array
  save_data::Bool
  AtomicSim() = new()
end

mutable struct Exchange <: MicroEnergy
   A::Array{Float64, 1}
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct Vector_Exchange <: MicroEnergy
   A::Array{Float64, 1}
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct ExchangeRKKY <: MicroEnergy
   sigma::Float64
   Delta::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct HeisenbergExchange <: Interaction
   Js::Array{Float64, 1}   #The length of Js should be equal to the number of neighbours
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct BulkDMI
   Dx::Float64
   Dy::Float64
   Dz::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct InterfacialDMI
   D::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end


mutable struct Anisotropy
   Ku::Array{Float64, 1}
   axis::Tuple
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct KagomeAnisotropy
   Ku::Float64
   ax1::Tuple
   ax2::Tuple
   ax3::Tuple
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct CubicAnisotropy
   Kc::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct Zeeman
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct TimeZeeman
   time_fun::Function
   init_field::Array{Float64, 1}
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end
