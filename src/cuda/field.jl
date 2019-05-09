using CUDAnative, CuArrays

function effective_field(zeeman::ZeemanGPU, sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  nxyz = sim.nxyz
  volume = sim.mesh.volume
  blocks_n, threads_n = sim.blocks, sim.threads
  @cuda blocks=blocks_n threads=threads_n zeeman_kernel!(spin, sim.field, zeeman.field_gpu, sim.energy, sim.Ms, volume, nxyz)
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

function effective_field(sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  fill!(sim.driver.field, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.driver.field .+= sim.field
  end
  return 0
end

function compute_system_energy(sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  sim.total_energy = 0
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
	interaction.total_energy = sum(sim.energy)
	sim.total_energy += interaction.total_energy
  end
  return 0
end


function compute_fields_to_gpu(sim::MicroSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
	copyto!(interaction.field, sim.field)
	copyto!(interaction.energy, sim.energy)
  end
  return 0
end
