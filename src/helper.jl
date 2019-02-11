function init_vector!(v::Array{Float64}, mesh::Mesh, init::Function)
  nxyz = mesh.nx*mesh.ny*mesh.nz
  dx,dy,dz = mesh.dx, mesh.dy, mesh.dz
  b = reshape(v, 3, nxyz)
	for i = 1:mesh.nx
		for j = 1:mesh.ny
			for k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        b[:, id] .=  init(i,j,k,dx*mesh.unit_length,dy*mesh.unit_length,dz*mesh.unit_length)
      end
    end
  end
end

function init_vector!(v::Array{Float64}, mesh::Mesh, init::Tuple{Number,Number,Number})
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
	b[1, :] .= init[1]
  b[2, :] .= init[2]
  b[3, :] .= init[3]
end

function init_vector!(v::Array{Float64}, mesh::Mesh, init::Array{Number})
	nxyz = mesh.nx*mesh.ny*mesh.nz
	b = reshape(v, 3, nxyz)
  if length(init) == 3
		b[1, :] .= init[1]
    b[2, :] .= init[2]
    b[3, :] .= init[3]
  elseif length(init) == length(v)
    v[:] .= init[:]
  else
    error("Init array length must be equal to 3*nx*ny*nz or 3")
  end
end

function normalise(a::Array, N::Int64)
  for i = 0:N-1
      j = 3*i+1
      length = sqrt(a[j]*a[j]+a[j+1]*a[j+1]+a[j+2]*a[j+2])
      if length > 0
        length = 1.0/length
				a[j] *= length
				a[j+1] *= length
				a[j+2] *= length
      end
    end
end

function omega_to_spin(omega::Array, spin::Array, spin_next::Array, N::Int64)
  #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
	#where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
  for i = 0:N-1
	  j = 3*i + 1
		w1 = omega[j]*0.5
		w2 = omega[j+1]*0.5
		w3 = omega[j+2]*0.5
		m1 = spin[j]
		m2 = spin[j+1]
		m3 = spin[j+2]
		r = 1 + w1*w1 + w2*w2 + w3*w3
		a11 = 1 + w1*w1 - w2*w2 - w3*w3
		a12 = 2*(w1*w2 - w3)
		a13 = 2*(w2 + w1*w3)
		a21 = 2*(w1*w2 + w3)
		a22 = 1 - w1*w1 + w2*w2 - w3*w3
		a23 = -2*(w1-w2*w3)
		a31 = 2*(-w2+w1*w3)
		a32 = 2*(w1+w2*w3)
		a33 = 1 - w1*w1 - w2*w2 + w3*w3
		spin_next[j] = (a11*m1 + a12*m2 + a13*m3)/r
		spin_next[j+1] = (a21*m1 + a22*m2 + a23*m3)/r
		spin_next[j+2] = (a31*m1 + a32*m2 + a33*m3)/r
  end
end

function compute_dmdt(m1::Array{Float64}, m2::Array{Float64}, N::Int64, dt::Float64)
  max_dmdt = 0.0
  for i = 0:N-1
    j = 3*i + 1
    mx = m1[j] - m2[j]
    my = m1[j+1] - m2[j+1]
    mz = m1[j+2] - m2[j+2]
    dmdt = sqrt(mx*mx + my*my + mz*mz)/dt
    if dmdt > max_dmdt
      max_dmdt = dmdt
    end
    return max_dmdt
  end
end
