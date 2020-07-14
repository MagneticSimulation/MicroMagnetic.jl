using CUDA

function effective_field(zeeman::ZeemanGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  volume = sim.mesh.volume
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n zeeman_kernel!(spin, sim.field, zeeman.cufield, sim.energy, sim.Ms, volume, nxyz)
  return nothing
end

function effective_field(stochastic::StochasticFieldGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
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
  k_B = 1.3806505e-23
  factor = 2*alpha*k_B/(mu_0*volume*gamma*dt)

  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n stochastic_field_kernel!(spin, sim.field, stochastic.eta, stochastic.T, sim.energy, sim.Ms, factor, volume, nxyz)
  return nothing
end

function effective_field(zeeman::TimeZeemanGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  volume = sim.mesh.volume
  blocks_n, threads_n = sim.blocks, sim.threads
  tx, ty, tz = zeeman.time_fun(t)
  @cuda blocks=blocks_n threads=threads_n time_zeeman_kernel!(spin, sim.field, zeeman.init_field, sim.energy, sim.Ms, volume, T(tx), T(ty), T(tz), nxyz)
  return nothing
end

function effective_field(exch::ExchangeGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_kernel!(spin, sim.field, sim.energy, sim.Ms, exch.A,
                    mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz, mesh.volume,
                    mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end

function effective_field(exch::Vector_ExchangeGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_vector_kernel!(spin, sim.field, sim.energy, sim.Ms, exch.A,
                    mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz, mesh.volume,
                    mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end

function effective_field(exch::ExchangeRKKYGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_kernel_rkky!(spin, sim.field, sim.energy, sim.Ms, exch.sigma/exch.Delta,
                    mesh.nx, mesh.ny, mesh.nz)
  return nothing
end

function effective_field(dmi::BulkDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n bulkdmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.Dx, dmi.Dy, dmi.Dz, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end

function effective_field(dmi::InterfacialDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n interfacial_dmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.D, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic)
  return nothing
end

function effective_field(dmi::SpatialInterfacialDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n spatial_interfacial_dmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.D, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic)
  return nothing
end


function effective_field(dmi::SpatialBulkDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n spatial_bulkdmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.D, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end


function effective_field(anis::AnisotropyGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
    blocks_n, threads_n = sim.blocks, sim.threads
    axis = anis.axis
    volume = sim.mesh.volume
    @cuda blocks=blocks_n threads=threads_n anisotropy_kernel!(spin, sim.field, sim.energy, anis.Ku,
                                         T(axis[1]), T(axis[2]),T(axis[3]),
                                         sim.Ms, volume, sim.nxyz)
    return nothing
end

function effective_field(anis::CubicAnisotropyGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
    blocks_n, threads_n = sim.blocks, sim.threads
    volume = sim.mesh.volume
    @cuda blocks=blocks_n threads=threads_n cubic_anisotropy_kernel!(spin, sim.field, sim.energy, anis.Kc,
                                         sim.Ms, volume, sim.nxyz)
    return nothing
end

function effective_field(anis::TitledCubicAnisotropyGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  blocks_n, threads_n = sim.blocks, sim.threads
  volume = sim.mesh.volume
  axis = anis.axis
  @cuda blocks=blocks_n threads=threads_n titled_cubic_anisotropy_kernel!(spin, sim.field, sim.energy, anis.Kc,
                                       axis[1],axis[2],axis[3],axis[4],axis[5],axis[6],axis[7],axis[8],axis[9],
                                       sim.Ms, volume, sim.nxyz)
  return nothing
end

function effective_field(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  fill!(sim.driver.field, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.driver.field .+= sim.field
  end
  return 0
end

function compute_system_energy(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  sim.total_energy = 0
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
	interaction.total_energy = sum(sim.energy)
	sim.total_energy += interaction.total_energy
  end
  return 0
end


function compute_fields_to_gpu(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
	copyto!(interaction.field, sim.field)
	copyto!(interaction.energy, sim.energy)
  end
  return 0
end
