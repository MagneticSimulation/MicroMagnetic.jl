
function compute_exch(spin::Array, field::Array, energy::Array, ngbs::Array, J::Float64)
  fd = [0.0,0.0,0.0]
  for i=1:size(spin)[2]
    fd[:] = 0.0
    for j=1:6
      id = ngbs[j,i]
      if id>0
        fd[:]+=J*spin[:,id]
      end
    end
    field[:,i] += fd[:]
    energy[i] = -0.5 * dot(fd[:], spin[:,i])
  end
  return 0
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
  n = sim.nxyz
  H = [sim.Hx, sim.Hy, sim.Hz]
  b = reshape(sim.field, 3, sim.nxyz)
  for i=1:n
   b[:, i] = H[:]
  end
  if sim.J!=0.0
    #compute_anis(spin, sim.field, sim.energy, sim.ngbs, sim.J)
  end
  return 0
end
