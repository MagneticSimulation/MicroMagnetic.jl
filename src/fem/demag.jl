using Printf
using LinearAlgebra

mutable struct DemagFE{T<:AbstractFloat}
    D::AbstractSparseMatrix
    G::AbstractSparseMatrix
    K1::Any
    K2::Any
    B::AbstractMatrix{T}
    g1::AbstractArray{T,1}
    g2::AbstractArray{T,1}
    phi1::AbstractArray{T,1}
    phi2::AbstractArray{T,1}
    u1_bnd::AbstractArray{T,1}
    u2_bnd::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    H::Any
    name::String
    method::String
    using_hmatrix::Bool
    DemagFE{T}() where {T<:AbstractFloat} = new()
end


#compute the solid angle 
function solid_angle_single(p::Array{Float64,1}, x1::Array{Float64,1}, x2::Array{Float64,1},
                            x3::Array{Float64,1})
    x = x1 .- p
    y = x2 .- p
    z = x3 .- p

    d = det3(x, y, z)
    a, b, c = norm2(x), norm2(y), norm2(z)

    div = a * b * c + dot3(x, y) * c + dot3(x, z) * b + dot3(y, z) * a
    omega = atan(abs(d), div)

    if omega < 0
        omega += pi
    end

    return 2 * omega
end

function boundary_element(xp::Array{Float64,1}, x1::Array{Float64,1}, x2::Array{Float64,1},
                          x3::Array{Float64,1})
    res = zeros(3)
    omega = solid_angle_single(xp, x1, x2, x3)
    if omega == 0
        return res
    end

    rhov1 = x1 .- xp
    rhov2 = x2 .- xp
    rhov3 = x3 .- xp
    sv1 = x2 .- x1
    sv2 = x3 .- x2
    sv3 = x1 .- x3

    zetav = cross3(sv1, sv2)

    xi1 = sv1 ./ norm2(sv1)
    xi2 = sv2 ./ norm2(sv2)
    xi3 = sv3 ./ norm2(sv3)
    zetav = zetav ./ norm2(zetav)

    zeta = dot3(zetav, rhov1)
    etav1 = cross3(zetav, xi1)
    etav2 = cross3(zetav, xi2)
    etav3 = cross3(zetav, xi3)

    eta1 = dot3(etav1, rhov1)
    eta2 = dot3(etav2, rhov2)
    eta3 = dot3(etav3, rhov3)

    gamma1 = [dot3(xi2, xi1), dot3(xi2, xi2), dot3(xi2, xi3)]
    gamma2 = [dot3(xi3, xi1), dot3(xi3, xi2), dot3(xi3, xi3)]
    gamma3 = [dot3(xi1, xi1), dot3(xi1, xi2), dot3(xi1, xi3)]

    s1 = norm2(sv1)
    s2 = norm2(sv2)
    s3 = norm2(sv3)
    s = (s1 + s2 + s3) / 2.0
    area = sqrt(s * (s - s1) * (s - s2) * (s - s3))

    rho1 = norm2(rhov1)
    rho2 = norm2(rhov2)
    rho3 = norm2(rhov3)

    p = zeros(3)
    p[1] = log((rho1 + rho2 + s1) / (rho1 + rho2 - s1))
    p[2] = log((rho2 + rho3 + s2) / (rho2 + rho3 - s2))
    p[3] = log((rho3 + rho1 + s3) / (rho3 + rho1 - s3))

    if zeta < 0
        omega = -omega
    end

    res[1] = (eta2 * omega - zeta * dot3(gamma1, p)) * s2 / (8.0 * pi * area)
    res[2] = (eta3 * omega - zeta * dot3(gamma2, p)) * s3 / (8.0 * pi * area)
    res[3] = (eta1 * omega - zeta * dot3(gamma3, p)) * s1 / (8.0 * pi * area)

    return res
end

function compute_degrees(K::SparseMatrixCSC)
    m, n = size(K)
    row_nnz = zeros(Int, m)
    
    for col in 1:n
        for i in nzrange(K, col)
            row = rowvals(K)[i]
            row_nnz[row] += 1
        end
    end
    
    return row_nnz
end

#assemble the matrix D, G and K1, and K2
function assemble_matirx_DGK1K2(demag::DemagFE, sim::MicroSimFE; using_constraint=true)
    mesh = sim.mesh
    D = spzeros(mesh.number_nodes, 3 * mesh.number_nodes)
    G = spzeros(3 * mesh.number_nodes, mesh.number_nodes)
    K1 = spzeros(mesh.number_nodes, mesh.number_nodes)
    K2 = spzeros(mesh.number_nodes, mesh.number_nodes)
    Ms = sim.mu0_Ms ./ mu_0
    map_g2b = Array(mesh.map_g2b)

    for c in 1:(mesh.number_cells)
        k1 = mesh.cell_verts[1, c]
        k2 = mesh.cell_verts[2, c]
        k3 = mesh.cell_verts[3, c]
        k4 = mesh.cell_verts[4, c]

        v = (k1, k2, k3, k4)

        bcd = compute_bcd12(mesh.coordinates[:, k1], mesh.coordinates[:, k2],
                            mesh.coordinates[:, k3], mesh.coordinates[:, k4])

        for alpha in 1:4, beta in 1:4
            i = 3 * (alpha - 1)
            j = 3 * (beta - 1)

            fa = mesh.volumes[c] * (bcd[i + 1] * bcd[j + 1] +
                                    bcd[i + 2] * bcd[j + 2] +
                                    bcd[i + 3] * bcd[j + 3])

            K1[v[alpha], v[beta]] += fa

            if map_g2b[v[alpha]] < 0
                K2[v[alpha], v[beta]] += fa
            end

            for k in 1:3
                dg1 = bcd[j + k] * mesh.volumes[c] / 4.0
                dg2 = dg1 * Ms[c]

                D[v[beta], 3 * (v[alpha] - 1) + k] += dg2

                G[3 * (v[alpha] - 1) + k, v[beta]] += dg1
            end
        end
    end

    map_b2g = Array(mesh.map_b2g)
    for i in 1:(mesh.number_nodes_bnd)
        j = map_b2g[i]
        K2[j, j] = 1.0
    end

    if using_constraint
        regions = maximum(mesh.region_ids)
        N_rank = rank(K1)
        if size(K1, 1) - N_rank != regions
            @warn("The number of mesh regions are not equal to the number of connected region!!!")
        end

        degrees = compute_degrees(K1)

        N_constraints = min(size(K1, 1) - N_rank, regions)
        for id = 1:N_constraints
            mask = mesh.region_ids .== id
            selected = mesh.cell_verts[1, mask]
            I = argmax(degrees[selected])
            J = selected[I]
            
            K1[J, :] .= 0
            K1[:, J] .= 0
            K1[J, J] = 1

            D[J, :] .= 0 #make sure demag.g1[I] = 0
        end
        
    end

    if default_backend[] != CPU()
        demag.D = GPUSparseMatrixCSC[](D)
        demag.G = GPUSparseMatrixCSC[](G)
        demag.K1 = GPUSparseMatrixCSR[](K1)
        demag.K2 = GPUSparseMatrixCSR[](K2)
    else
        demag.D = D
        demag.G = G
        demag.K1 = K1
        demag.K2 = K2
    end

end


function ComputeVertBSA(mesh::FEMesh)
    vertbsa = zeros(mesh.number_nodes)
    for c in 1:(mesh.number_cells), j in 0:3
        k1 = mesh.cell_verts[(j + 0) % 4 + 1, c]
        k2 = mesh.cell_verts[(j + 1) % 4 + 1, c]
        k3 = mesh.cell_verts[(j + 2) % 4 + 1, c]
        k4 = mesh.cell_verts[(j + 3) % 4 + 1, c]

        v1 = mesh.coordinates[:, k1]
        v2 = mesh.coordinates[:, k2]
        v3 = mesh.coordinates[:, k3]
        v4 = mesh.coordinates[:, k4]

        omega = solid_angle_single(v1, v2, v3, v4)

        vertbsa[k1] += omega
    end

    return vertbsa
end

function assemble_matrix_B(demag::DemagFE, sim::MicroSimFE)
    mesh = sim.mesh
    B = zeros(mesh.number_nodes_bnd, mesh.number_nodes_bnd)
    map_g2b = Array(mesh.map_g2b)

    Threads.@threads for n in 1:(mesh.number_nodes)

        #if the node is not at the boundary, do nothing
        if map_g2b[n] < 0
            continue
        end

        #local v = mesh.coordinates[:, n]

        for f in 1:(mesh.number_faces_bnd)
            local i = mesh.face_verts[1, f]
            local j = mesh.face_verts[2, f]
            local k = mesh.face_verts[3, f]

            # skip if the node belongs to the triangle
            if n == i || n == j || n == k
                continue
            end

            local be = boundary_element(mesh.coordinates[:, n], mesh.coordinates[:, i],
                                        mesh.coordinates[:, j], mesh.coordinates[:, k])

            B[map_g2b[n], map_g2b[i]] += be[1]
            B[map_g2b[n], map_g2b[j]] += be[2]
            B[map_g2b[n], map_g2b[k]] += be[3]
        end
    end

    vertbsa = ComputeVertBSA(mesh)

    for n in 1:(mesh.number_nodes)
        if map_g2b[n] > 0
            tmp = vertbsa[n] / (4.0 * pi) - 1.0
            B[map_g2b[n], map_g2b[n]] += tmp
        end
    end

    if default_backend[] != CPU()
        demag.B = GPUMatrix[](B)
    else
        demag.B = B
    end
    
end

function init_demag(sim::MicroSimFE, method; kwargs...)
    mesh = sim.mesh

    T = Float[]
    demag = DemagFE{T}()
    demag.method = method

    demag.g1 = create_zeros(mesh.number_nodes)
    demag.g2 = create_zeros(mesh.number_nodes)
    demag.phi1 = create_zeros(mesh.number_nodes)
    demag.phi2 = create_zeros(mesh.number_nodes)

    demag.u1_bnd = create_zeros(mesh.number_nodes_bnd)
    demag.u2_bnd = create_zeros(mesh.number_nodes_bnd)

    demag.field = create_zeros(3 * mesh.number_nodes)
    demag.energy = create_zeros(mesh.number_nodes)

    assemble_matirx_DGK1K2(demag, sim)

    try
        demag.K1 = cholesky(demag.K1) #reuse the factorization to speed up the calculation
    catch
        try
            demag.K1 = lu(demag.K1)
        catch
            demag.K1 = qr(demag.K1)
        end
    end

    try
        demag.K2 = lu(demag.K2)
    catch
        demag.K2 = qr(demag.K2)
    end

    demag.using_hmatrix = false

    if demag.method == "bem_hmatrix"
        demag.using_hmatrix = true
        @info("BEM HMatrix demag is added")
        #HMatrixGPU.set_backend("cpu")
        #demag.H = init_Hmatrix(sim, mesh; kwargs...)
        #@info HMatrixGPU.info(demag.H)
    elseif demag.method == "bem"
        @info "BEM demag is added"
        @info "Start to build BEM matrix, will take a while ..."
        time = @elapsed assemble_matrix_B(demag, sim)
        @info @sprintf("Building the BEM matrix is finished, it takes %0.2f seconds.", time)
    else
        @info("Supported methods for demag: bem_hmatrix and bem")
    end

    return demag
end