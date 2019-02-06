module SpinDynamics

include("mesh.jl")

type SimData
  spin::Array{Float64}
  field::Array{Float64}
  energy::Array{Float64}
  nxyz::Int64
  name::ASCIIString
  alpha::Float64
  J::Float64
  K::Float64
  Hx::Float64
  Hy::Float64
  Hz::Float64
end

function create_sim(mesh::Mesh; name="dyn")
  nxyz = mesh.nx*mesh.ny*mesh.nz
  spin = zeros(Float64,3,nxyz)
  field = zeros(Float64,3,nxyz)
  energy = zeros(Float64,nxyz)
  return SimData(spin, field, energy, nxyz, name, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

include("ode.jl")

end
