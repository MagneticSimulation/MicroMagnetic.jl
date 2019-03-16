struct Mesh
  dx::Float64
  dy::Float64
  dz::Float64
  nx::Int64
  ny::Int64
  nz::Int64
  unit_length::Float64
  volume::Float64
  ngbs::Array{Int64}
  pbc::String
end


mutable struct Dopri5
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::Array{Float64}
   omega_t::Array{Float64}
   dw_dt::Array{Float64}
   ks::Array{Float64}
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

mutable struct SimData
  mesh::Mesh
  ode::Dopri5
  saver::DataSaver
  spin::Array{Float64}
  prespin::Array{Float64}
  field::Array{Float64}
  energy::Array{Float64}
  Ms::Array{Float64}
  nxyz::Int64
  name::String
  alpha::Float64
  gamma::Float64
  precession::Bool
  interactions::Array
end

mutable struct Exchange
   A::Float64
   field::Array{Float64}
   energy::Array{Float64}
   name::String
end

mutable struct BulkDMI
   D::Float64
   field::Array{Float64}
   energy::Array{Float64}
   name::String
end

mutable struct Anisotropy
   Ku::Array{Float64}
   axis::Tuple
   field::Array{Float64}
   energy::Array{Float64}
   name::String
end

mutable struct Zeeman
   Hx::Float64
   Hy::Float64
   Hz::Float64
   field::Array{Float64}
   energy::Array{Float64}
   name::String
end
