function effective_field(zee::Zeeman, sim::MicroSimFEM, spin::Array{Float64, 1}, t::Float64)
    mu0 = 4*pi*1e-7
    #nxyz = sim.nxyz
    field = zee.field
    #volume = sim.mesh.volume
    #for i = 1:nxyz
    #  j = 3*(i-1)
    #zee.energy[i] = -mu0*sim.Ms[i]*volume*(spin[j+1]*field[j+1] + spin[j+2]*field[j+2] + spin[j+3]*field[j+3])
    #end
  end