abstract type Driver end
abstract type Mesh end
abstract type AbstractSim end
abstract type MicroEnergy end
abstract type Interaction end
abstract type Integrator end
abstract type AbstractDopri5 <:Integrator end

struct FDMesh <: Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  volume::Float64
  ngbs::Array{Int64, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

struct CubicMesh <: Mesh
  a::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  ngbs::Array{Int64, 2}
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
end

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
end

mutable struct AtomicSim <: AbstractSim
  mesh::CubicMesh
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
end


mutable struct Exchange <: MicroEnergy
   A::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct HeisenbergExchange <: Interaction
   J::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end

mutable struct BulkDMI
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

mutable struct Zeeman
   Hx::Float64
   Hy::Float64
   Hz::Float64
   field::Array{Float64, 1}
   energy::Array{Float64, 1}
   name::String
end
