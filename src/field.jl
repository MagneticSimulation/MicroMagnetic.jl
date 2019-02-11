mutable struct Exchange
   A::Float64
	 field::Array{Float64}
   energy::Array{Float64}
   name::String
end

mutable struct BulkDMI
   D::Float64
	 field::Array{Float64}
   energy::Array{Float64}
   name::String
end

mutable struct Anisotropy
   Ku::Array{Float64}
   axis::Tuple
	 field::Array{Float64}
   energy::Array{Float64}
   name::String
end

mutable struct Zeeman
   Hx::Float64
   Hy::Float64
   Hz::Float64
	 field::Array{Float64}
   energy::Array{Float64}
   name::String
end

function effective_field(zee::Zeeman, sim::SimData, spin::Array{Float64}, t::Float64)
  nxyz = sim.nxyz
	b = reshape(zee.field, 3, nxyz)
	for i = 1:nxyz
	  b[1, i] = zee.Hx
    b[2, i] = zee.Hy
    b[3, i] = zee.Hz
	end
end

function effective_field(anis::Anisotropy, sim::SimData, spin::Array{Float64}, t::Float64)
	mu0 = 4.0*pi*1e-7
  mesh = sim.mesh
  nxyz = sim.nxyz
	field = anis.field
  energy = anis.energy
  Ms = sim.Ms
  Ku = anis.Ku
  axis = anis.axis
	for i = 1:nxyz
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
		field[k+1] = 2*Ku[i]*spin[k+1]*axis[1]*Ms_inv
		field[k+2] = 2*Ku[i]*spin[k+2]*axis[2]*Ms_inv
		field[k+3] = 2*Ku[i]*spin[k+3]*axis[3]*Ms_inv
		energy[i] = Ku[i]*(1.0-sa*sa)
	end

end

function effective_field(exch::Exchange, sim::SimData, spin::Array{Float64}, t::Float64)
  mu0 = 4.0*pi*1e-7
  mesh = sim.mesh
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
        fx += A[j]*(spin[k+1]-spin[3*i-2])
        fy += A[j]*(spin[k+2]-spin[3*i-1])
        fz += A[j]*(spin[k+3]-spin[3*i])
      end
    end
		Ms_inv = 1.0/(Ms[i]*mu0)
    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])
    field[3*i-2] = fx*Ms_inv
    field[3*i-1] = fy*Ms_inv
    field[3*i] = fz*Ms_inv
  end
end

function effective_field(dmi::BulkDMI, sim::SimData, spin::Array{Float64}, t::Float64)
  mu0 = 4*pi*1e-7
  mesh = sim.mesh
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

    Ms_inv = 1.0/(Ms[i]*mu0)
    energy[i] = -0.5*(fx*spin[3*i-2] + fy*spin[3*i-1] + fz*spin[3*i])
    field[3*i-2] = fx*Ms_inv
    field[3*i-1] = fy*Ms_inv
    field[3*i] = fz*Ms_inv
  end
end


function effective_field(sim::SimData, spin::Array, t::Float64)
  fill!(sim.field, 0.0)
  for interaction in sim.interactions
    effective_field(interaction, sim, spin, t)
    sim.field[:] += interaction.field[:]
  end
  return 0
end
