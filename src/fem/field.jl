function assemble_anis_matirx(anis::AnisotropyFE, sim::MicroSimFE)
    mesh = sim.mesh
    Ku = anis.Ku
    K_mat = anis.K_matrix
    axis = anis.axis
    mu0_Ms = sim.mu0_Ms

    for c in 1:(mesh.number_cells)
        k1 = mesh.cell_verts[1, c]
        k2 = mesh.cell_verts[2, c]
        k3 = mesh.cell_verts[3, c]
        k4 = mesh.cell_verts[4, c]

        v = ((k1, k2, k3, k4) .- 1) .* 3

        for alpha in 1:4, beta in 1:4
            coef = -2 * Ku[c] * mesh.volumes[c] * (1 + kronecker_delta(alpha, beta)) /
                   (20 * mu0_Ms[c])

            for i in 1:3, j in 1:3
                fk = coef * axis[3 * (c - 1) + i] * axis[3 * (c - 1) + j]

                K_mat[v[alpha] + i, v[beta] + j] += fk
            end
        end
    end

    if default_backend[] != CPU()
        anis.K_matrix = GPUSparseMatrixCSC[](K_mat)
    end
    
end

"""
  assemble_exch_matirx(exch::ExchangeFE, sim::MicroSimFE)
"""
function assemble_exch_matirx(exch::ExchangeFE, sim::MicroSimFE)
    mesh = sim.mesh
    unit_length = mesh.unit_length
    K = exch.K_matrix
    A = exch.A
    mu0_Ms = sim.mu0_Ms

    for c in 1:(mesh.number_cells)
        k1 = mesh.cell_verts[1, c]
        k2 = mesh.cell_verts[2, c]
        k3 = mesh.cell_verts[3, c]
        k4 = mesh.cell_verts[4, c]

        v = ((k1, k2, k3, k4) .- 1) .* 3

        bcd = compute_bcd12(mesh.coordinates[:, k1], mesh.coordinates[:, k2],
                            mesh.coordinates[:, k3], mesh.coordinates[:, k4])

        coef = 2 * (A[c] / unit_length) * mesh.volumes[c] / mu0_Ms[c] / unit_length

        for alpha in 0:3, beta in 0:3
            i = 3 * alpha + 1
            j = 3 * beta + 1

            fa = coef *
                 (bcd[i] * bcd[j] + bcd[i + 1] * bcd[j + 1] + bcd[i + 2] * bcd[j + 2])

            for k in 1:3
                K[v[alpha + 1] + k, v[beta + 1] + k] += fa
            end
        end
    end
    if default_backend[] != CPU()
        exch.K_matrix = GPUSparseMatrixCSC[](K)
    end
end



function effective_field(zee::Zeeman, sim::MicroSimFE, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    v_coeff = sim.mesh.unit_length^3

    back = default_backend[]
    zeeman_fe_kernel!(back, groupsize[])(spin, zee.field, zee.energy, sim.L_mu, T(v_coeff);
                                      ndrange=N)
    return nothing
end


function effective_field(anis::Union{AnisotropyFE, ExchangeFE}, sim::MicroSimFE,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    mesh = sim.mesh
    N = sim.n_total

    mul!(anis.field, anis.K_matrix, spin)
    anis.field .*= mesh.L_inv_neg

    v_coeff = 0.5 * mesh.unit_length^3
    back = default_backend[]
    zeeman_fe_kernel!(back, groupsize[])(spin, anis.field, anis.energy, sim.L_mu, T(v_coeff);
                                      ndrange=N)

    return nothing
end


function effective_field(demag::DemagFE, sim::MicroSimFE, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}

    mesh = sim.mesh

    # step1: phi1 = D*m
    mul!(demag.g1, demag.D, spin)

    # step2: solve K1*phi1 = g1
    demag.phi1 = demag.K1 \ demag.g1

    # step3: u1_bnd = U1*phi1
    #mul!(demag.u1_bnd, demag.U1, demag.phi1)
    demag.u1_bnd .= demag.phi1[mesh.map_b2g]

    # step4: u2_bnd = B*u1_bnd
    if demag.using_hmatrix
        #HMatrixGPU.mul!(demag.u2_bnd, demag.H, demag.u1_bnd)
    else
        mul!(demag.u2_bnd, demag.B, demag.u1_bnd)
    end

    # step5: g2 = U2*u2_bnd
    #mul!(demag.g2, demag.U2, demag.u2_bnd)
    demag.g2[mesh.map_b2g] .= demag.u2_bnd

    # step6: solve K2*phi2 = g2
    demag.phi2 = demag.K2 \ demag.g2

    # step7: phi = phi1 + phi2
    demag.phi1 .+= demag.phi2

    # step8: F = G*phi
    mul!(demag.field, demag.G, demag.phi1)

    demag.field .*= mesh.L_inv_neg

    v_coeff = 0.5 * mesh.unit_length^3
    backend = default_backend[]
    zeeman_fe_kernel!(backend, groupsize[])(spin, demag.field, demag.energy, sim.L_mu, T(v_coeff); ndrange=sim.n_total)

end