mutable struct Exchange
   A::Float64
	 field::Array{Float64}
   energy::Array{Float64}
end

mutable struct BulkDMI
   D::Float64
	 field::Array{Float64}
   energy::Array{Float64}
end

mutable struct Zeeman
   Hx::Float64
   Hy::Float64
   Hz::Float64
	 field::Array{Float64}
   energy::Array{Float64}
   name::String
end

function effective_field(zee::Zeeman, sim::SimData, mesh::Mesh)
  nxyz = sim.nxyz
	b = reshape(zee.field, 3, nxyz)
	for i = 1:nxyz
	  b[1, i] = zee.Hx
    b[2, i] = zee.Hy
    b[3, i] = zee.Hz
	end

end

function effective_field(exch::Exchange, sim::SimData, mesh::Mesh)
  mu0 = 4*pi*1e-7
  dx = mesh.dx*mesh.unit_length
  dy = mesh.dy*mesh.unit_length
  dz = mesh.dz*mesh.unit_length
  ngbs = mesh.ngbs
  nxyz = sim.nxyz
  field = exch.field
  energy = exch.energy
  Ms = sim.Ms
	ax = 2.0 * exch.A / (dx * dx)
  ay = 2.0 * exch.A / (dy * dy)
	az = 2.0 * exch.A / (dz * dz)
  A = (ax, ax, ay, ay, az, az)
  for i = 1:nxyz
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
    for j=1:6
      id = ngbs[j,i]
      if id>0 && Ms[id]>0
        k = 3*(id-1)
        fx += A[j]*spin[k+1]
        fy += A[j]*spin[k+2]
        fz += A[j]*spin[k+3]
      end
    end
    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])
    field[3*i-2] = fx/(mu0*Ms[i])
    field[3*i-1] = fy/(mu0*Ms[i])
    field[3*i] = fz/(mu0*Ms[i])
  end
end

function effective_field(dmi::BulkDMI, sim::SimData, mesh::Mesh)
  mu0 = 4*pi*1e-7
  dx = mesh.dx*mesh.unit_length
  dy = mesh.dy*mesh.unit_length
  dz = mesh.dz*mesh.unit_length
  ngbs = mesh.ngbs
  nxyz = sim.nxyz
  field = dmi.field
  energy = dmi.energy
  Ms = sim.Ms
  D = sim.D
  Ds = (D/dx, D/dx, D/dy, D/dy, D/dz, D/dz)
  ax = (1.0,-1.0, 0.0, 0.0, 0.0, 0.0)
  ay = (0.0, 0.0, 1.0,-1.0, 0.0, 0.0)
  az = (0.0, 0.0, 0.0, 0.0, 1.0,-1.0)

  for i = 1:nxyz
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
      id = ngbs[1,j]
      if id>0 && Ms[id]>0
        k = 3*(id-1)+1
			  fx += Ds[j]*cross_x(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
			  fy += Ds[j]*cross_y(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
			  fz += Ds[j]*cross_z(ax[j],ay[j],az[j],spin[k],spin[k+1],spin[k+2]);
      end
    end

    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])
    field[3*i-2] = fx/(mu0*Ms[i])
    field[3*i-1] = fy/(mu0*Ms[i])
    field[3*i] = fz/(mu0*Ms[i])
  end
end


function compute_anis(spin::Array, field::Array, energy::Array, Ku::Array, axis::Array)
  for i=1:size(spin)[2]
    sa = spin[1,i]*axis[1]+spin[2,i]*axis[2]+spin[3,i]*axis[3]
    field[1,i] = 2*Ku[i]*axis[1]
    field[2,i] = 2*Ku[i]*axis[2]
    field[3,i] = 2*Ku[i]*axis[3]
    energy[i] = -Ku[i]*sa*sa
  end
  return 0
end


function effective_field(sim::SimData, spin::Array, t::Float64)
  fill!(sim.field, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, sim.mesh)
    sim.field[:] += interaction.field[:]
  end
  return 0
end
