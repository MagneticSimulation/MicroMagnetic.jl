struct EnergyMinimization <: Driver
  m::Array{Float64, 1}
  m_next::Array{Float64, 1}
  tau::Float64
  nsteps::Int64
  nfevals::Int64
end

struct LLG <: Driver
  m::Array{Float64, 1}
  m_next::Array{Float64, 1}
  tau::Float64
  nsteps::Int64
  nfevals::Int64
end


function create_driver(driver::String, nxyz::Int64) #TODO: FIX ME
    if driver=="SDM"
        m = zeros(Float64,3*nxyz)
        m_next = zeros(Float64,3*nxyz)
        return EnergyMinimization(m, m_next, 0.0, 0, 0)
    end
    return nothing
end
