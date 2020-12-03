const mu_0 = 4 * pi * 1e-7
const mu_B = 9.27400949e-24
const k_B = 1.3806505e-23
const c_e = 1.602176565e-19
const eV = 1.602176565e-19
const meV = 1.602176565e-22
const m_e = 9.10938291e-31
const g_e = 2.0023193043737
const h_bar = 1.05457172647e-34
const gamma = g_e * mu_B / h_bar
const mu_s_1 = g_e * mu_B * 1.0  # for S=1, 1.856952823077189e-23
const h_bar_gamma = h_bar * gamma
const mT = 0.001 / (4*pi*1e-7)


abstract type AbstractSim end
abstract type AbstractSimGPU <:AbstractSim end
abstract type MonteCarloEnergy end

mutable struct ZeemanMC{TF<:AbstractFloat}
    H::CuArray{TF, 1}#h1,h2,h3
end

mutable struct Anis6Fold2DMC{TF<:AbstractFloat} <: MonteCarloEnergy
    K::CuArray{TF, 1}#k1,k2
end

mutable struct ExchangeMC{TF<:AbstractFloat} <: MonteCarloEnergy
    J::CuArray{TF, 1} #nearest exchange constant
end

mutable struct MonteCarlo{TF<:AbstractFloat} <:AbstractSimGPU
  mesh::Mesh
  exch::MonteCarloEnergy
  # dmi::MonteCarloEnergy
  zee::ZeemanMC
  anis::MonteCarloEnergy
  # saver::DataSaver
  spin::CuArray{TF, 1}
  nextspin::CuArray{TF, 1}
  rnd::CuArray{TF, 1}
  energy::CuArray{TF, 1}
  energy2::CuArray{TF, 1}
  delta_E::CuArray{TF, 1}
  EnMean::Float64
  En2Mean::Float64
  nxyz::Int64
  steps::Int64
  name::String
  T::Float64
  Qdens::CuArray{TF, 1}
  spin2::CuArray{TF, 1}
  Q::Float64
  mxMean::TF
  myMean::TF
  mzMean::TF
  mx2Mean::TF
  my2Mean::TF
  mz2Mean::TF
  MonteCarlo{T}() where {T<:AbstractFloat} = new()
end

