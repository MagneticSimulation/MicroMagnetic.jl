using CuArrays

struct FDMeshGPU <: Mesh
  dx::FloatGPU
  dy::FloatGPU
  dz::FloatGPU
  nx::Int64
  ny::Int64
  nz::Int64
  nxyz::Int64
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
  volume::FloatGPU
end

function FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc=(false, false, false))
  nxyz = nx*ny*nz
  volume = dx*dy*dz
  return FDMeshGPU(FloatGPU(dx), FloatGPU(dy), FloatGPU(dz), nx, ny, nz, nxyz, pbc[1], pbc[2], pbc[3], FloatGPU(volume))
end

mutable struct Dopri5GPU
   tol::Float64
   t::Float64
   step::Float64
   step_next::Float64
   facmax::Float64
   facmin::Float64
   safety::Float64
   nsteps::Int64
   nfevals::Int64
   omega::Array{FloatGPU, 1}
   omega_t::Array{FloatGPU, 1}
   dw_dt::Array{FloatGPU, 1}
   ks::Array{FloatGPU, 2}
   rhs_fun::Function
   succeed::Bool
end

mutable struct MicroSimGPU <: AbstractSim
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
