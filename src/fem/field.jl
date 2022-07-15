
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


function assemble_exch_matirx(exch::ExchangeFEM, sim::MicroSimFEM)
  
  mesh = sim.mesh
  unit_length = mesh.unit_length
  K = exch.K_matrix
  A = exch.A
  Ms = sim.Ms


  for c = 1 : mesh.number_cells

    k1 = mesh.cell_verts[1, c]
    k2 = mesh.cell_verts[2, c]
    k3 = mesh.cell_verts[3, c]
    k4 = mesh.cell_verts[4, c]

    v = ((k1, k2, k3, k4) .- 1) .* 3
    #println(v,   K.m,  K.n)

    bcd = compute_bcd12(mesh.coordinates[:, k1], 
                        mesh.coordinates[:, k2],
                        mesh.coordinates[:, k3],
                        mesh.coordinates[:, k4])

    coef = 2 * (A[c] / unit_length) * mesh.volumes[c] / (mu_0 * Ms[c]) / unit_length;

    for alpha = 0:3, beta = 0:3
        i = 3 * alpha + 1
        j = 3 * beta + 1

        fa = coef*(bcd[i] * bcd[j] + bcd[i + 1] * bcd[j + 1] + bcd[i + 2] * bcd[j + 2]);

        for k = 1:3
          K[v[alpha+1] + k, v[beta+1] + k] += fa
        end
  
    end
  end

end

function effective_field(exch::ExchangeFEM, sim::MicroSimFEM, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7

  mesh = sim.mesh

  if !exch.K_matrix_initialized
    exch.K_matrix_initialized = true
    assemble_exch_matirx(exch, sim)
  end
  
  field = exch.K_matrix*spin

  exch.field .= (field .* mesh.L_inv_neg)
  #volume = sim.mesh.volume
  #for i = 1:nxyz
  #  j = 3*(i-1)
  #zee.energy[i] = -mu0*sim.Ms[i]*volume*(spin[j+1]*field[j+1] + spin[j+2]*field[j+2] + spin[j+3]*field[j+3])
  #end
end


