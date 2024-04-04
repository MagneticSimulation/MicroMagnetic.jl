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

function effective_field(exch::ExchangeGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_kernel!(spin, sim.field, sim.energy, sim.Ms, exch.A,
                    mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz, mesh.volume,
                    mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
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

function effective_field(exch::Vector_ExchangeGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n exchange_vector_kernel!(spin, sim.field, sim.energy, sim.Ms, exch.A,
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

function effective_field(dmi::BulkDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n bulkdmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.Dx, dmi.Dy, dmi.Dz, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
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

function effective_field(dmi::SpatialInterfacialDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n spatial_interfacial_dmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.D, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic)
  return nothing
end


function effective_field(dmi::SpatialBulkDMIGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  n_total = sim.n_total
  mesh = sim.mesh
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n spatial_bulkdmi_kernel!(spin, sim.field, sim.energy, sim.Ms,
                    dmi.D, mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny,
                    mesh.nz, mesh.volume, mesh.xperiodic, mesh.yperiodic, mesh.zperiodic)
  return nothing
end


"""
    effective_field(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

GPU version effective_field function.

For example:

```julia
#sim = Sim(GPUmesh)
effective_field(sim, sim.spin, 0.0)
```

After running this function, the effective field is calculated and saved in sim.driver.field, which is a CPU array.
"""
function effective_field(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  fill!(sim.driver.field, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.driver.field .+= sim.field
  end
  return 0
end

"""
    compute_system_energy(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

Calculate the total energy of sim, which will be saved in sim.total_energy.

Parameters is same with function effective_field above.
"""
function compute_system_energy(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  sim.total_energy = 0
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
	interaction.total_energy = sum(sim.energy)
	sim.total_energy += interaction.total_energy
  end
  return 0
end

"""
    compute_fields_to_gpu(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

Calculate the effective field for every interactions of sim.
Fields will be saved in sim.interaction.field, and energy will be saved in sim.interaction.energy, which both are CPU array.
"""
function compute_fields_to_gpu(sim::AbstractSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
	copyto!(interaction.field, sim.field)
	copyto!(interaction.energy, sim.energy)
  end
  return 0
end
