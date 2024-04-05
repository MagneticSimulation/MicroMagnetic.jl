using CUDA

function effective_field(stochastic::StochasticFieldGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  volume = sim.mesh.volume
  integrator = sim.driver.ode

  if integrator.nsteps > stochastic.nsteps
      randn!(stochastic.eta)
      stochastic.nsteps = integrator.nsteps
  end

  mu0 = 4*pi*1e-7
  dt = integrator.step
  gamma = sim.driver.gamma
  alpha = sim.driver.alpha
  k_B = stochastic.k_B
  factor = 2*alpha*k_B/(mu_0*volume*gamma*dt)

  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n stochastic_field_kernel!(spin, sim.field, stochastic.eta, stochastic.T, sim.energy, sim.Ms, factor, volume, n_total)
  return nothing
end

function effective_field(exch::ExchangeAnistropyGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_anistropy_kernel!(spin, sim.field, sim.energy, sim.Ms, exch.kea,
                    mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz, mesh.volume,
                    mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end

function effective_field(exch::ExchangeRKKYGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_kernel_rkky!(spin, sim.field, sim.energy, sim.Ms, exch.J/mesh.dz, mesh.volume,
                    mesh.nx, mesh.ny, mesh.nz)
  return nothing
end


function effective_field(dmi::InterlayerDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n interlayer_dmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.Dx, dmi.Dy, dmi.Dz, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end

function effective_field(dmi::InterfacialDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n interfacial_dmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.D, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic)
  return nothing
end
