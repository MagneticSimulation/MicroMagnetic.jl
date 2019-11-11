using Random

abstract type MCInteraction end

mutable struct ExchangeMC <: MCInteraction
   Jx::Float64 #nearest exchange constant
   Jy::Float64
   Jz::Float64
   Jx1::Float64 #next-nearest exchange constant
   Jy1::Float64
   Jz1::Float64
   name::String
end

mutable struct DMI_MC <: MCInteraction
   Dx::Float64 #nearest DMI constant
   Dy::Float64
   Dz::Float64
   Dx1::Float64 #next-nearest DMI constant
   Dy1::Float64
   Dz1::Float64
   name::String
end

mutable struct Anisotropy_MC <: MCInteraction
   Ku::Float64 #nearest DMI constant
   axis_ux::Float64
   axis_uy::Float64
   axis_uz::Float64
   Kc::Float64
   axis_cx::Float64
   axis_cy::Float64
   axis_cz::Float64
   name::String
end


mutable struct MonteCarloNew{TF<:AbstractFloat} <:AbstractSimGPU
  mesh::Mesh
  mu_s::Float64
  interactions::Array{Any, 1}
  saver::DataSaver
  spin::CuArray{TF, 1}
  nextspin::CuArray{TF, 1}
  rnd::CuArray{TF, 1}
  energy::CuArray{TF, 1}
  delta_E::CuArray{TF, 1}
  total_energy::TF
  nxyz::Int64
  steps::Int64
  name::String
  T::Float64
  MonteCarloNew{T}() where {T<:AbstractFloat} = new()
end

function MonteCarloNew(mesh::Mesh; name="mc")
    Float = _cuda_using_double.x ? Float64 : Float32
    sim = MonteCarloNew{Float}()
    sim.mesh = mesh
    nxyz = mesh.nx*mesh.ny*mesh.nz
    sim.nxyz = nxyz

    sim.spin = CuArrays.zeros(Float, 3*nxyz)
    sim.nextspin = CuArrays.zeros(Float,3*nxyz)
    sim.rnd = CuArrays.zeros(Float,3*nxyz)
    sim.energy = CuArrays.zeros(Float,nxyz)
    sim.delta_E = CuArrays.zeros(Float,nxyz)
    sim.steps = 0
    sim.name = name
    sim.T = 300

    headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
    units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
    results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> o.total_energy, average_m]
    saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
    sim.saver = saver

    return sim
end

function init_m0(sim::MonteCarloNew, m0::Any; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end
  copyto!(sim.spin, spin)
  return true
end

function run_single_step_triangular(sim::MonteCarlo, bias::Int64)
  F = _cuda_using_double.x ? Float64 : Float32
  blk, thr = CuArrays.cudims(sim.nxyz)
  @cuda blocks=blk threads=thr run_step_triangular_kernel!(sim.spin, sim.nextspin,
                               sim.rnd, sim.energy, sim.mesh.ngbs,
                               F(sim.J), F(sim.J1), F(sim.D), F(sim.D1), sim.bulk_dmi, F(sim.Ku),
                               F(sim.Hx), F(sim.Hy), F(sim.Hz), F(sim.T),
                               sim.mesh.nx, sim.mesh.ny, bias)
  return nothing
end

function run_single_step_cubic(sim::MonteCarlo, bias::Int64)
  F = _cuda_using_double.x ? Float64 : Float32
  blk, thr = CuArrays.cudims(sim.nxyz)
  mesh = sim.mesh
  @cuda blocks=blk threads=thr run_step_cubic_kernel!(sim.spin, sim.nextspin,
                               sim.rnd, sim.energy, mesh.ngbs, mesh.nngbs,
                               F(sim.J), F(sim.J1), F(sim.D), F(sim.D1), sim.bulk_dmi,
                               F(sim.Ku), F(sim.Kc),
                               F(sim.Hx), F(sim.Hy), F(sim.Hz), F(sim.T),
                               mesh.nx, mesh.ny, mesh.nz, bias)
  return nothing
end


function run_step_triangular(sim::MonteCarlo)

  uniform_random_sphere(sim.nextspin, sim.rnd, sim.nxyz)
  run_single_step_triangular(sim, 0)
  run_single_step_triangular(sim, 1)
  run_single_step_triangular(sim, 2)
  sim.steps += 1

  return  nothing
end

function run_step_cubic(sim::MonteCarlo)

  uniform_random_sphere(sim.nextspin, sim.rnd, sim.nxyz)
  run_single_step_cubic(sim, 0)
  run_single_step_cubic(sim, 1)
  run_single_step_cubic(sim, 2)
  sim.steps += 1

  return  nothing
end

function compute_triangular_energy(sim::MonteCarlo)
    F = _cuda_using_double.x ? Float64 : Float32
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr total_energy_triangular_kernel!(sim.spin, sim.energy,
                                 sim.mesh.ngbs,
                                 F(sim.J), F(sim.J1), F(sim.D), F(sim.D1),
                                 sim.bulk_dmi, F(sim.Ku),
                                 F(sim.Hx), F(sim.Hy), F(sim.Hz),
                                 sim.mesh.nxyz)
    return sum(sim.energy)*k_B
end

function compute_cubic_energy(sim::MonteCarlo)
    F = _cuda_using_double.x ? Float64 : Float32
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr total_energy_cubic_kernel!(sim.spin, sim.energy,
                                 sim.mesh.ngbs, sim.mesh.nngbs,
                                 F(sim.J), F(sim.J1), F(sim.D), F(sim.D1),
                                 sim.bulk_dmi, F(sim.Ku),
                                 F(sim.Hx), F(sim.Hy), F(sim.Hz),
                                 sim.mesh.nxyz)
    return sum(sim.energy)*k_B
end
