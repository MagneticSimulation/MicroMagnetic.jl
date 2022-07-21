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

function assemble_anis_matirx(anis::AnisotropyFEM, sim::MicroSimFEM)
  
  mesh = sim.mesh
  Ku = anis.Ku
  K_mat = anis.K_matrix
  axis = anis.axis
  Ms = sim.Ms

  for c = 1 : mesh.number_cells

    k1 = mesh.cell_verts[1, c]
    k2 = mesh.cell_verts[2, c]
    k3 = mesh.cell_verts[3, c]
    k4 = mesh.cell_verts[4, c]

    v = ((k1, k2, k3, k4) .- 1) .* 3

    for alpha = 1:4, beta = 1:4

        coef = -2 * Ku[c] * mesh.volumes[c]*(1 + kronecker_delta(alpha, beta)) / (20 * mu_0 * Ms[c])

        for i = 1:3, j = 1:3
          
          fk = coef*axis[3*(c-1)+i]*axis[3*(c-1)+j] 

          K_mat[v[alpha]+i, v[beta]+j] += fk
        end
  
    end
  end

end

function effective_field(zee::Zeeman, sim::MicroSimFEM, spin::Array{Float64, 1}, t::Float64)
  
  mu0 = 4*pi*1e-7
  field = zee.field
  mesh = sim.mesh

  f_Ms = field .* sim.L_mu

  v_coeff = mesh.unit_length*mesh.unit_length*mesh.unit_length
  
  for i = 1:mesh.number_nodes
    j = 3*(i-1)
    zee.energy[i] = -mu0*(spin[j+1]*f_Ms[j+1] + spin[j+2]*f_Ms[j+2] + spin[j+3]*f_Ms[j+3])*v_coeff
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

  f_Ms = exch.field .* sim.L_mu

  v_coeff = mesh.unit_length*mesh.unit_length*mesh.unit_length
  
  for i = 1:sim.n_nodes
    j = 3*(i-1)
    exch.energy[i] = -0.5*mu0*(spin[j+1]*f_Ms[j+1] + spin[j+2]*f_Ms[j+2] + spin[j+3]*f_Ms[j+3])*v_coeff
  end
end


function effective_field(anis::AnisotropyFEM, sim::MicroSimFEM, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7

  mesh = sim.mesh

  if !anis.K_matrix_initialized
    anis.K_matrix_initialized = true
    assemble_anis_matirx(anis, sim)
  end
  
  field = anis.K_matrix*spin

  anis.field .= (field .* mesh.L_inv_neg)

  f_Ms = anis.field .* sim.L_mu

  v_coeff = mesh.unit_length*mesh.unit_length*mesh.unit_length
  
  for i = 1:sim.n_nodes
    j = 3*(i-1)
    anis.energy[i] = -0.5*mu0*(spin[j+1]*f_Ms[j+1] + spin[j+2]*f_Ms[j+2] + spin[j+3]*f_Ms[j+3])*v_coeff
  end
end


function effective_field(demag::DemagFEM, sim::MicroSimFEM, spin::Array{Float64, 1}, t::Float64)
  mu0 = 4*pi*1e-7
  mesh = sim.mesh

  if !demag.matrices_assembled
    demag.matrices_assembled = true
    assemble_matirx_DGK1K2(demag, sim)
    assmeble_matrix_B(demag, sim)
  end
  
  # step1: phi1 = D*m
  mul!(demag.g1,demag.D,spin)

  # step2: solve K1*phi1 = g1
  demag.phi1 = demag.K1 \ demag.g1

  # step3: u1_bnd = U1*phi1
  #mul!(demag.u1_bnd, demag.U1, demag.phi1)
  for i = 1:mesh.number_nodes_bnd
    id = mesh.map_b2g[i]
    demag.u1_bnd[i] = demag.phi1[id]
  end
  
  # step4: u2_bnd = B*u1_bnd
  mul!(demag.u2_bnd, demag.B, demag.u1_bnd)

  # step5: g2 = U2*u2_bnd
  #mul!(demag.g2, demag.U2, demag.u2_bnd)
  for i = 1:mesh.number_nodes_bnd
    id = mesh.map_b2g[i]
    demag.g2[id] = demag.u2_bnd[i]
  end
  
  # step6: solve K2*phi2 = g2
  demag.phi2 = demag.K2 \ demag.g2

  # step7: phi = phi1 + phi2
  demag.phi1 .+= demag.phi2

  # step8: F = G*phi
  mul!(demag.field, demag.G, demag.phi1)

  demag.field .*= mesh.L_inv_neg

  f_Ms = demag.field .* sim.L_mu

  v_coeff = mesh.unit_length*mesh.unit_length*mesh.unit_length
  
  for i = 1:sim.n_nodes
    j = 3*(i-1)
    demag.energy[i] = -0.5*mu0*(spin[j+1]*f_Ms[j+1] + spin[j+2]*f_Ms[j+2] + spin[j+3]*f_Ms[j+3])*v_coeff
  end
end

