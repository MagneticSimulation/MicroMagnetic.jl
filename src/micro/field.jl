function effective_field(zee::Zeeman, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7
  n_total = sim.n_total
  field = zee.field
  volume = sim.mesh.volume
  for i = 1:n_total
    j = 3*(i-1)
    zee.energy[i] = -mu0*sim.Ms[i]*volume*(spin[j+1]*field[j+1] + spin[j+2]*field[j+2] + spin[j+3]*field[j+3])
  end
end

function effective_field(zee::TimeZeeman, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7
  n_total = sim.n_total
  field = zee.field
  b = reshape(field, 3, n_total)
  b0 = reshape(zee.init_field, 3, n_total)
  volume = sim.mesh.volume
  tx, ty, tz = zee.time_fun(t)
  b[1, :] = b0[1, :] * tx
  b[2, :] = b0[2, :] * ty
  b[3, :] = b0[3, :] * tz
  for i = 1:n_total
    j = 3*(i-1)
    zee.energy[i] = -mu0*sim.Ms[i]*volume*(spin[j+1]*field[j+1] + spin[j+2]*field[j+2] + spin[j+3]*field[j+3])
  end
end

function effective_field(anis::Anisotropy, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4.0*pi*1e-7
  mesh = sim.mesh
  n_total = sim.n_total
  field = anis.field
  energy = anis.energy
  Ms = sim.Ms
  Ku = anis.Ku
  axis = anis.axis
  for i = 1:n_total
    if Ms[i] == 0.0
      energy[i] = 0.0
      field[3*i-2] = 0.0
      field[3*i-1] = 0.0
      field[3*i] = 0.0
      continue
    end
    k = 3*(i-1)
    sa = spin[k+1]*axis[1]+spin[k+2]*axis[2]+spin[k+3]*axis[3]
    Ms_inv = 1.0/(Ms[i]*mu0)
    field[k+1] = 2*Ku[i]*sa*axis[1]*Ms_inv
    field[k+2] = 2*Ku[i]*sa*axis[2]*Ms_inv
    field[k+3] = 2*Ku[i]*sa*axis[3]*Ms_inv
    energy[i] = Ku[i]*(1.0-sa*sa)*mesh.volume
  end

end

function effective_field(anis::CubicAnisotropy, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4.0*pi*1e-7
  mesh = sim.mesh
  n_total = sim.n_total
  field = anis.field
  energy = anis.energy
  Ms = sim.Ms
  Kc = anis.Kc
  axis = anis.axis
  for i = 1:n_total
    if Ms[i] == 0.0
      energy[i] = 0.0
      field[3*i-2] = 0.0
      field[3*i-1] = 0.0
      field[3*i] = 0.0
      continue
    end
    j = 3*(i-1)
    Ms_inv = 4*Kc/(Ms[i]*mu0)

    mxp = axis[1]*spin[j+1] + axis[2]*spin[j+2] + axis[3]*spin[j+3]
    myp = axis[4]*spin[j+1] + axis[5]*spin[j+2] + axis[6]*spin[j+3]
    mzp = axis[7]*spin[j+1] + axis[8]*spin[j+2] + axis[9]*spin[j+3]
    mxp3 = mxp*mxp*mxp
    myp3 = myp*myp*myp
    mzp3 = mzp*mzp*mzp
    field[j+1] = Ms_inv*(mxp3*axis[1]+myp3*axis[4]+mzp3*axis[7])
    field[j+2] = Ms_inv*(mxp3*axis[2]+myp3*axis[5]+mzp3*axis[8])
    field[j+3] = Ms_inv*(mxp3*axis[3]+myp3*axis[6]+mzp3*axis[9])
    energy[i] = -Kc*(mxp*mxp3+myp*myp3+mzp*mzp3)*mesh.volume
  end
end

function effective_field(exch::Exchange, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4.0*pi*1e-7
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
  ax = 2.0  / (dx * dx)
  ay = 2.0  / (dy * dy)
  az = 2.0  / (dz * dz)
  nabla = (ax, ax, ay, ay, az, az)

  for index = 1:n_total
    i = 3*index - 2
    if Ms[index] == 0.0
        energy[index] = 0.0
        field[i] = 0.0
        field[i+1] = 0.0
        field[i+2] = 0.0
      continue
    end
    fx, fy, fz = 0.0, 0.0, 0.0
    for j=1:6
      id = ngbs[j,index]
      if id>0 && Ms[id]>0
        k = 3*id-2
        fx += A[index]*nabla[j]*(spin[k]-spin[i])
        fy += A[index]*nabla[j]*(spin[k+1]-spin[i+1])
        fz += A[index]*nabla[j]*(spin[k+2]-spin[i+2])
      end
    end
    Ms_inv = 1.0/(Ms[index]*mu0)
    energy[index] = -0.5*(fx*spin[i] + fy*spin[i+1] + fz*spin[i+2])*mesh.volume
    field[i] = fx*Ms_inv
    field[i+1] = fy*Ms_inv
    field[i+2] = fz*Ms_inv
  end
end

function effective_field(exch::Vector_Exchange, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4.0*pi*1e-7
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
  ax = 2.0  / (dx * dx)
  ay = 2.0  / (dy * dy)
  az = 2.0  / (dz * dz)
  nabla = (ax, ax, ay, ay, az, az)

  for index = 1:n_total
    i = 3*index - 2
    if Ms[index] == 0.0
        energy[index] = 0.0
        field[i] = 0.0
        field[i+1] = 0.0
        field[i+2] = 0.0
        continue
    end

    fx, fy, fz = 0.0, 0.0, 0.0
    for j=1:2
      id = ngbs[j,index]
      if id>0 && Ms[id]>0
        k = 3*id-2
        fx += A[i]*nabla[j]*(spin[k]-spin[i])
        fy += A[i]*nabla[j]*(spin[k+1]-spin[i+1])
        fz += A[i]*nabla[j]*(spin[k+2]-spin[i+2])
      end
    end
    for j=3:4
      id = ngbs[j,index]
      if id>0 && Ms[id]>0
        k = 3*id-2
        fx += A[i+1]*nabla[j]*(spin[k]-spin[i])
        fy += A[i+1]*nabla[j]*(spin[k+1]-spin[i+1])
        fz += A[i+1]*nabla[j]*(spin[k+2]-spin[i+2])
      end
    end
    for j=5:6
      id = ngbs[j,index]
      if id>0 && Ms[id]>0
        k = 3*id-2
        fx += A[i+2]*nabla[j]*(spin[k]-spin[i])
        fy += A[i+2]*nabla[j]*(spin[k+1]-spin[i+1])
        fz += A[i+2]*nabla[j]*(spin[k+2]-spin[i+2])
      end
    end
    Ms_inv = 1.0/(Ms[index]*mu0)
    energy[index] = -0.5*(fx*spin[i] + fy*spin[i+1] + fz*spin[i+2])*mesh.volume
    field[i] = fx*Ms_inv
    field[i+1] = fy*Ms_inv
    field[i+2] = fz*Ms_inv
  end
end


function effective_field(exch::ExchangeRKKY, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4.0*pi*1e-7
  mesh = sim.mesh
  field = exch.field
  energy = exch.energy
  Ms = sim.Ms
  sigma = exch.J/mesh.dz
  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  for i=1:nx, j=1:ny
      id1 = (j-1) * nx + i
      id2 = (nz-1) * nx*ny + (j-1) * nx + i
      k1 = 3*id1-2
      k2 = 3*id2-2
      mbx = spin[k1]
      mby = spin[k1+1]
      mbz = spin[k1+2]

      mtx = spin[k2]
      mty = spin[k2+1]
      mtz = spin[k2+2]

      if Ms[id1] > 0 && Ms[id2] > 0
          Ms_inv = 1.0/(Ms[id1]*mu0)
          field[k1] = sigma*Ms_inv*mtx
          field[k1+1] = sigma*Ms_inv*mty
          field[k1+2] = sigma*Ms_inv*mtz
          energy[id1] = 0.5*sigma*(1-mtx*mbx-mty*mby-mtz*mbz)

          Ms_inv = 1.0/(Ms[id2]*mu0)
          field[k2] = sigma*Ms_inv*mbx
          field[k2+1] = sigma*Ms_inv*mby
          field[k2+2] = sigma*Ms_inv*mbz
          energy[id2] = energy[id1]
      end

  end

end

function effective_field(dmi::BulkDMI, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7
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
  Ds = (Dx/dx, Dx/dx, Dy/dy, Dy/dy, Dz/dz, Dz/dz)
  ax = (1.0,-1.0, 0.0, 0.0, 0.0, 0.0)
  ay = (0.0, 0.0, 1.0,-1.0, 0.0, 0.0)
  az = (0.0, 0.0, 0.0, 0.0, 1.0,-1.0)

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
      id = ngbs[j,i]
      if id>0 && Ms[id]>0
        k = 3*(id-1)+1
        fx += Ds[j]*cross_x(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
        fy += Ds[j]*cross_y(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
        fz += Ds[j]*cross_z(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
      end
    end

    Ms_inv = 1.0/(Ms[i]*mu0)
    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])*mesh.volume
    field[3*i-2] = fx*Ms_inv
    field[3*i-1] = fy*Ms_inv
    field[3*i] = fz*Ms_inv
  end
end


function effective_field(dmi::DMI, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7
  mesh = sim.mesh
  dx = mesh.dx
  dy = mesh.dy
  dz = mesh.dz
  ngbs = mesh.ngbs
  n_total = sim.n_total
  field = dmi.field
  energy = dmi.energy
  Ms = sim.Ms
  a = -dmi.gamma_A
  b = -0.5*(dmi.alpha_S+dmi.alpha_A-dmi.beta_S+dmi.beta_A)
  c = -0.5*(-dmi.alpha_S+dmi.alpha_A+dmi.beta_S+dmi.beta_A)

  for i = 1:n_total

    if Ms[i] == 0.0
      energy[i] = 0.0
      field[3*i-2] = 0.0
      field[3*i-1] = 0.0
      field[3*i] = 0.0
      continue
    end

    #x-direction
    i1 = ngbs[1,i]
    i2 = ngbs[2,i]
    factor = i1*i2>0 ? 1/(2*dx) : 1/dx
    i1 < 0 && (i1 = i)
    i2 < 0 && (i2 = i)
    j1 = 3*i1-2
    j2 = 3*i2-2
    PxSx = (spin[j2] - spin[j1]) * factor
    PxSy = (spin[j2+1] - spin[j1+1]) * factor
    PxSz = (spin[j2+2] - spin[j1+2]) * factor

    #y-direction
    i1 = ngbs[3,i]
    i2 = ngbs[4,i]
    factor = i1*i2>0 ? 1/(2*dy) : 1/dy
    i1 < 0 && (i1 = i)
    i2 < 0 && (i2 = i)
    j1 = 3*i1-2
    j2 = 3*i2-2
    PySx = (spin[j2] - spin[j1]) * factor
    PySy = (spin[j2+1] - spin[j1+1]) * factor
    PySz = (spin[j2+2] - spin[j1+2]) * factor


    #z-direction
    i1 = ngbs[5,i]
    i2 = ngbs[6,i]
    factor = i1*i2>0 ? 1/(2*dz) : 1/dz
    i1 < 0 && (i1 = i)
    i2 < 0 && (i2 = i)
    j1 = 3*i1-2
    j2 = 3*i2-2
    PzSx = (spin[j2] - spin[j1]) * factor
    PzSy = (spin[j2+1] - spin[j1+1]) * factor
    PzSz = (spin[j2+2] - spin[j1+2]) * factor

    fx = -dmi.xi_12*PxSy - dmi.xi_13*PxSz + dmi.xi_21*PySy + 2*c*PySz - 2*a*PzSy + dmi.xi_31*PzSz
    fy = dmi.xi_12*PxSx - 2*b*PxSz - dmi.xi_21*PySx - dmi.xi_23*PySz + 2*a*PzSx + dmi.xi_32*PzSz
    fz = dmi.xi_13*PxSx + 2*b*PxSy - 2*c*PySx + dmi.xi_23*PySy - dmi.xi_31*PzSx - dmi.xi_32*PzSy

    Ms_inv = 1.0/(Ms[i]*mu0)
    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])*mesh.volume
    field[3*i-2] = fx*Ms_inv
    field[3*i-1] = fy*Ms_inv
    field[3*i] = fz*Ms_inv
  end
end


function effective_field(dmi::InterfacialDMI, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7
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
  Ds = (D/dx, D/dx, D/dy, D/dy)
  ax = (0.0, 0.0, -1.0, 1.0) #Dij = D r_ij x z
  ay = (1.0,-1.0, 0.0, 0.0)
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
      id = ngbs[j,i]
      if id>0 && Ms[id]>0
        k = 3*(id-1)+1
        fx += Ds[j]*cross_x(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
        fy += Ds[j]*cross_y(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
        fz += Ds[j]*cross_z(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
      end
    end

    Ms_inv = 1.0/(Ms[i]*mu0)
    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])*mesh.volume
    field[3*i-2] = fx*Ms_inv
    field[3*i-1] = fy*Ms_inv
    field[3*i] = fz*Ms_inv
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
function effective_field(sim::AbstractSim, spin::Array{Float64, 1}, t::Float64)
  fill!(sim.field, 0.0)
  fill!(sim.energy, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.field .+= interaction.field
    sim.energy .+= interaction.energy
  end
  return 0
end

function compute_system_energy(sim::AbstractSim, spin::Array{Float64, 1}, t::Float64)
  #sim.total_energy = 0
  fill!(sim.energy, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.energy .+= interaction.energy
  end
  return 0
end
