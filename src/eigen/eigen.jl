#using SparseArrays

export build_matrix

function build_matrix(sim; gamma=2.21e5, sparse=false, alpha=0.01)
    @info("Building matrix ...")
    spin = zeros(AbstractFloat, 3 * sim.n_total)
    m0 = Array(sim.spin)

    # add perturbation to static magnetization m0
    for i in 1:(sim.n_total)
        x = 3 * (i - 1) + 1
        R = rotation_matrix(m0[x], m0[x + 1], m0[x + 2])
        v = ones(AbstractFloat, 3)
        v[1] = Epsilon(2 * i - 1, 1.0)
        v[2] = Epsilon(2 * i, 1.0)
        spin[x] = R[1,1]*v[1] + R[1, 2]*v[2] + R[1, 3]*v[3]
        spin[x+1] = R[2,1]*v[1] + R[2, 2]*v[2] + R[2, 3]*v[3]
        spin[x+2] = R[3,1]*v[1] + R[3, 2]*v[2] + R[3, 3]*v[3]
    end

    for i in 1:length(spin)
        spin[i] = simplify(spin[i])
    end

    @info("  Calculating the effective field ...")
    # compute the effective field H_eff
    MicroMagnetic.effective_field(sim, spin)

    H_eff = sim.field
    for i in 1:length(H_eff)
        H_eff[i] = simplify(H_eff[i])
    end

    @info("  Calculating m x H ...")
    dm_dt = zeros(AbstractFloat, 3 * sim.n_total)
    for i in 1:(sim.n_total)
        x = 3 * (i - 1) + 1
        dm_dt[x] = -gamma*cross_x(spin[x], spin[x+1], spin[x+2], H_eff[x], H_eff[x+1], H_eff[x+2])
        dm_dt[x+1] = -gamma*cross_y(spin[x], spin[x+1], spin[x+2], H_eff[x], H_eff[x+1], H_eff[x+2])
        dm_dt[x+2] = -gamma*cross_z(spin[x], spin[x+1], spin[x+2], H_eff[x], H_eff[x+1], H_eff[x+2])
    end

    for i in 1:length(dm_dt)
        dm_dt[i] = simplify(dm_dt[i])
    end

    
    # we use spin to store the local dm_dt
    for i in 1:(sim.n_total)
        x = 3 * (i - 1) + 1
        R = rotation_matrix_inverse(m0[x], m0[x + 1], m0[x + 2])

        v1 = dm_dt[x]
        v2 = dm_dt[x+1]
        v3 = dm_dt[x+2]

        spin[x] = R[1,1]*v1 + R[1, 2]*v2 + R[1, 3]*v3
        spin[x+1] = R[2,1]*v1 + R[2, 2]*v2 + R[2, 3]*v3
        spin[x+2] = R[3,1]*v1 + R[3, 2]*v2 + R[3, 3]*v3
    end

    N = sim.n_total
    #matrix = sparse ? spzeros(2N, 2N) : zeros(2N, 2N)
    matrix = zeros(2N, 2N)
    @info("  Collecting terms ...")
    for i=1:N
        x = 3 * (i - 1) + 1
        terms = collect_terms(spin[x])
        for j=1:2N
            v = get(terms, j, 0)
            if abs(v) > eps()
                matrix[2*i-1, j] = v
            end
        end

        terms = collect_terms(spin[x+1])
        for j=1:2N
            v = get(terms, j, 0)
            if abs(v) > eps()
                matrix[2*i, j] = v
            end
        end
    end

    return matrix
end