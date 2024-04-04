function effective_field(zee::Zeeman, sim::MicroSim, spin, t::Float64)
  n_total = sim.n_total
  volume = sim.mesh.volume

  groupsize = 512
  kernel! = zeeman_kernel!(backend[], groupsize)
  kernel!(spin, zee.field, zee.energy, sim.Ms, volume, ndrange=n_total)

  KernelAbstractions.synchronize(backend[])
  return nothing
end

function effective_field(zee::TimeZeeman, sim::MicroSim, spin, t::Float64)
  n_total = sim.n_total
  volume = sim.mesh.volume
  tx, ty, tz = zee.time_fun(t)
  T = single_precision.x ? Float32 : Float64
  groupsize = 512
  kernel! = time_zeeman_kernel!(backend[], groupsize)
  kernel!(spin, zee.field, zee.init_field, zee.energy, sim.Ms, volume, T(tx), T(ty), T(tz), ndrange=n_total)
  KernelAbstractions.synchronize(backend[])
  return nothing
end

function effective_field(anis::Anisotropy, sim::MicroSim, spin, t::Float64)
  n_total = sim.n_total
  volume = sim.mesh.volume

  T = single_precision.x ? Float32 : Float64

  groupsize = 512
  kernel! = anisotropy_kernel!(backend[], groupsize)
  kernel!(spin, anis.field, anis.energy, anis.Ku, T(anis.axis[1]), T(anis.axis[2]), T(anis.axis[3]), sim.Ms, volume, ndrange=n_total)

  KernelAbstractions.synchronize(backend[])
  return nothing
end

function effective_field(anis::CubicAnisotropy, sim::MicroSim, spin, t::Float64)
  n_total = sim.n_total
  volume = sim.mesh.volume

  T = single_precision.x ? Float32 : Float64

  groupsize = 512
  kernel! = cubic_anisotropy_kernel!(backend[], groupsize)
  a1x, a1y, a1z = anis.axis1
  a2x, a2y, a2z = anis.axis2
  a3x, a3y, a3z = anis.axis3
  kernel!(spin, anis.field, anis.energy, anis.Kc, T(a1x), T(a1y), T(a1z), T(a2x), T(a2y), T(a2z),
    T(a3x), T(a3y), T(a3z), sim.Ms, volume, ndrange=n_total)

  KernelAbstractions.synchronize(backend[])
  return nothing
end

function effective_field(exch::Exchange, sim::MicroSim, spin::Array{Float64,1}, t::Float64)
  mu0 = 4.0 * pi * 1e-7
  mesh = sim.mesh
  dx = mesh.dx
  dy = mesh.dy
  dz = mesh.dz
  ngbs = mesh.ngbs
  n_total = sim.n_total
  field = exch.field
  energy = exch.energy
  A = exch.A
  Ms = sim.Ms
  ax = 2.0 / (dx * dx)
  ay = 2.0 / (dy * dy)
  az = 2.0 / (dz * dz)
  nabla = (ax, ax, ay, ay, az, az)

  for index = 1:n_total
    i = 3 * index - 2
    if Ms[index] == 0.0
      energy[index] = 0.0
      field[i] = 0.0
      field[i+1] = 0.0
      field[i+2] = 0.0
      continue
    end
    fx, fy, fz = 0.0, 0.0, 0.0
    for j = 1:6
      id = ngbs[j, index]
      if id > 0 && Ms[id] > 0
        k = 3 * id - 2
        fx += A[index] * nabla[j] * (spin[k] - spin[i])
        fy += A[index] * nabla[j] * (spin[k+1] - spin[i+1])
        fz += A[index] * nabla[j] * (spin[k+2] - spin[i+2])
      end
    end
    Ms_inv = 1.0 / (Ms[index] * mu0)
    energy[index] = -0.5 * (fx * spin[i] + fy * spin[i+1] + fz * spin[i+2]) * mesh.volume
    field[i] = fx * Ms_inv
    field[i+1] = fy * Ms_inv
    field[i+2] = fz * Ms_inv
  end
end

function effective_field(exch::VectorExchange, sim::MicroSim, spin::Array{Float64,1}, t::Float64)
  mu0 = 4.0 * pi * 1e-7
  mesh = sim.mesh
  dx = mesh.dx
  dy = mesh.dy
  dz = mesh.dz
  ngbs = mesh.ngbs
  n_total = sim.n_total
  field = exch.field
  energy = exch.energy
  A = exch.A
  Ms = sim.Ms
  ax = 2.0 / (dx * dx)
  ay = 2.0 / (dy * dy)
  az = 2.0 / (dz * dz)
  nabla = (ax, ax, ay, ay, az, az)

  for index = 1:n_total
    i = 3 * index - 2
    if Ms[index] == 0.0
      energy[index] = 0.0
      field[i] = 0.0
      field[i+1] = 0.0
      field[i+2] = 0.0
      continue
    end

    fx, fy, fz = 0.0, 0.0, 0.0
    for j = 1:2
      id = ngbs[j, index]
      if id > 0 && Ms[id] > 0
        k = 3 * id - 2
        fx += A[i] * nabla[j] * (spin[k] - spin[i])
        fy += A[i] * nabla[j] * (spin[k+1] - spin[i+1])
        fz += A[i] * nabla[j] * (spin[k+2] - spin[i+2])
      end
    end
    for j = 3:4
      id = ngbs[j, index]
      if id > 0 && Ms[id] > 0
        k = 3 * id - 2
        fx += A[i+1] * nabla[j] * (spin[k] - spin[i])
        fy += A[i+1] * nabla[j] * (spin[k+1] - spin[i+1])
        fz += A[i+1] * nabla[j] * (spin[k+2] - spin[i+2])
      end
    end
    for j = 5:6
      id = ngbs[j, index]
      if id > 0 && Ms[id] > 0
        k = 3 * id - 2
        fx += A[i+2] * nabla[j] * (spin[k] - spin[i])
        fy += A[i+2] * nabla[j] * (spin[k+1] - spin[i+1])
        fz += A[i+2] * nabla[j] * (spin[k+2] - spin[i+2])
      end
    end
    Ms_inv = 1.0 / (Ms[index] * mu0)
    energy[index] = -0.5 * (fx * spin[i] + fy * spin[i+1] + fz * spin[i+2]) * mesh.volume
    field[i] = fx * Ms_inv
    field[i+1] = fy * Ms_inv
    field[i+2] = fz * Ms_inv
  end
end


function effective_field(exch::ExchangeRKKY, sim::MicroSim, spin::Array{Float64,1}, t::Float64)
  mu0 = 4.0 * pi * 1e-7
  mesh = sim.mesh
  field = exch.field
  energy = exch.energy
  Ms = sim.Ms
  sigma = exch.J / mesh.dz
  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  for i = 1:nx, j = 1:ny
    id1 = (j - 1) * nx + i
    id2 = (nz - 1) * nx * ny + (j - 1) * nx + i
    k1 = 3 * id1 - 2
    k2 = 3 * id2 - 2
    mbx = spin[k1]
    mby = spin[k1+1]
    mbz = spin[k1+2]

    mtx = spin[k2]
    mty = spin[k2+1]
    mtz = spin[k2+2]

    if Ms[id1] > 0 && Ms[id2] > 0
      Ms_inv = 1.0 / (Ms[id1] * mu0)
      field[k1] = sigma * Ms_inv * mtx
      field[k1+1] = sigma * Ms_inv * mty
      field[k1+2] = sigma * Ms_inv * mtz
      energy[id1] = 0.5 * sigma * (1 - mtx * mbx - mty * mby - mtz * mbz)

      Ms_inv = 1.0 / (Ms[id2] * mu0)
      field[k2] = sigma * Ms_inv * mbx
      field[k2+1] = sigma * Ms_inv * mby
      field[k2+2] = sigma * Ms_inv * mbz
      energy[id2] = energy[id1]
    end

  end

end

function effective_field(dmi::BulkDMI, sim::MicroSim, spin::Array{Float64,1}, t::Float64)
  mu0 = 4 * pi * 1e-7
  mesh = sim.mesh
  dx = mesh.dx
  dy = mesh.dy
  dz = mesh.dz
  ngbs = mesh.ngbs
  n_total = sim.n_total
  field = dmi.field
  energy = dmi.energy
  Ms = sim.Ms
  Dx, Dy, Dz = dmi.Dx, dmi.Dy, dmi.Dz
  Ds = (Dx / dx, Dx / dx, Dy / dy, Dy / dy, Dz / dz, Dz / dz)
  ax = (1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
  ay = (0.0, 0.0, 1.0, -1.0, 0.0, 0.0)
  az = (0.0, 0.0, 0.0, 0.0, 1.0, -1.0)

  Threads.@threads for i = 1:n_total
    if Ms[i] == 0.0
      energy[i] = 0.0
      field[3*i-2] = 0.0
      field[3*i-1] = 0.0
      field[3*i] = 0.0
      continue
    end
    fx = 0.0
    fy = 0.0
    fz = 0.0

    for j = 1:6
      id = ngbs[j, i]
      if id > 0 && Ms[id] > 0
        k = 3 * (id - 1) + 1
        fx += Ds[j] * cross_x(ax[j], ay[j], az[j], spin[k], spin[k+1], spin[k+2])
        fy += Ds[j] * cross_y(ax[j], ay[j], az[j], spin[k], spin[k+1], spin[k+2])
        fz += Ds[j] * cross_z(ax[j], ay[j], az[j], spin[k], spin[k+1], spin[k+2])
      end
    end

    Ms_inv = 1.0 / (Ms[i] * mu0)
    energy[i] = -0.5 * (fx * spin[3*i-2] + fy * spin[3*i-1] + fz * spin[3*i]) * mesh.volume
    field[3*i-2] = fx * Ms_inv
    field[3*i-1] = fy * Ms_inv
    field[3*i] = fz * Ms_inv
  end
end


function effective_field(dmi::InterfacialDMI, sim::MicroSim, spin::Array{Float64,1}, t::Float64)
  mu0 = 4 * pi * 1e-7
  mesh = sim.mesh
  dx = mesh.dx
  dy = mesh.dy
  dz = mesh.dz
  ngbs = mesh.ngbs
  n_total = sim.n_total
  field = dmi.field
  energy = dmi.energy
  Ms = sim.Ms
  D = dmi.D
  Ds = (D / dx, D / dx, D / dy, D / dy)
  ax = (0.0, 0.0, -1.0, 1.0) #Dij = D r_ij x z
  ay = (1.0, -1.0, 0.0, 0.0)
  az = (0.0, 0.0, 0.0, 0.0)

  Threads.@threads for i = 1:n_total
    if Ms[i] == 0.0
      energy[i] = 0.0
      field[3*i-2] = 0.0
      field[3*i-1] = 0.0
      field[3*i] = 0.0
      continue
    end
    fx = 0.0
    fy = 0.0
    fz = 0.0

    for j = 1:4
      id = ngbs[j, i]
      if id > 0 && Ms[id] > 0
        k = 3 * (id - 1) + 1
        fx += Ds[j] * cross_x(ax[j], ay[j], az[j], spin[k], spin[k+1], spin[k+2])
        fy += Ds[j] * cross_y(ax[j], ay[j], az[j], spin[k], spin[k+1], spin[k+2])
        fz += Ds[j] * cross_z(ax[j], ay[j], az[j], spin[k], spin[k+1], spin[k+2])
      end
    end

    Ms_inv = 1.0 / (Ms[i] * mu0)
    energy[i] = -0.5 * (fx * spin[3*i-2] + fy * spin[3*i-1] + fz * spin[3*i]) * mesh.volume
    field[3*i-2] = fx * Ms_inv
    field[3*i-1] = fy * Ms_inv
    field[3*i] = fz * Ms_inv
  end
end

"""
    effective_field(sim::AbstractSim, spin::Array{Float64, 1}, t::Float64)

Calculate the effective field of a CPU Sim.

Parameters:

  sim : AbstractSim struct whose field is to be calculated.

  spin : 1-d array that matches sim.mesh. 

  t : Time in second. This term is used when TimeZeeman is added into sim.

For example:

  ```julia
    #sim = Sim(CPUmesh)
    effective_field(sim, sim.spin, 0.0)
  ```

After running this function, the effective field is calculated and saved in sim.field, 
and the total energy is saved in sim.energy.
"""
function effective_field(sim::AbstractSim, spin, t::Float64)
  fill!(sim.field, 0.0)
  fill!(sim.energy, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t::Float64)
    sim.field .+= interaction.field
    sim.energy .+= interaction.energy
  end
  return 0
end

function compute_system_energy(sim::AbstractSim, spin::Array{Float64,1}, t::Float64)
  #sim.total_energy = 0
  fill!(sim.energy, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.energy .+= interaction.energy
  end
  return 0
end
