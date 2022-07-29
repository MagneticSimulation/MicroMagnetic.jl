mutable struct DemagFEM
    D::SparseMatrixCSC{}
    G::SparseMatrixCSC{}
    K1::Any
    K2::Any
    B::Matrix{Float64}
    g1::Array{Float64, 1}
    g2::Array{Float64, 1}
    phi1::Array{Float64, 1}
    phi2::Array{Float64, 1}
    u1_bnd::Array{Float64, 1}
    u2_bnd::Array{Float64, 1}
    field::Array{Float64, 1}
    energy::Array{Float64, 1}
    name::String
    matrices_assembled::Bool
    DemagFEM() = new()
  end

#compute the solid angle 
function solid_angle_single(p::Array{Float64,1}, x1::Array{Float64,1}, 
                            x2::Array{Float64,1}, x3::Array{Float64,1})
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

function boundary_element(xp::Array{Float64,1}, x1::Array{Float64,1}, 
                          x2::Array{Float64,1}, x3::Array{Float64,1})
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


    zeta = dot3(zetav, rhov1);
    etav1 = cross3(zetav, xi1);
    etav2 = cross3(zetav, xi2);
    etav3 = cross3(zetav, xi3);

    eta1 = dot3(etav1, rhov1);
    eta2 = dot3(etav2, rhov2);
    eta3 = dot3(etav3, rhov3);

    gamma1 = [dot3(xi2, xi1), dot3(xi2, xi2), dot3(xi2, xi3)]
    gamma2 = [dot3(xi3, xi1), dot3(xi3, xi2), dot3(xi3, xi3)]
    gamma3 = [dot3(xi1, xi1), dot3(xi1, xi2), dot3(xi1, xi3)]

    s1 = norm2(sv1);
    s2 = norm2(sv2);
    s3 = norm2(sv3);
    s = (s1 + s2 + s3) / 2.0;
    area = sqrt(s * (s - s1) * (s - s2) * (s - s3));

    rho1 = norm2(rhov1);
    rho2 = norm2(rhov2);
    rho3 = norm2(rhov3);

    p = zeros(3)
    p[1] = log((rho1 + rho2 + s1) / (rho1 + rho2 - s1));
    p[2] = log((rho2 + rho3 + s2) / (rho2 + rho3 - s2));
    p[3] = log((rho3 + rho1 + s3) / (rho3 + rho1 - s3));

    if zeta < 0
        omega = -omega;
    end

    res[1] = (eta2 * omega - zeta * dot3(gamma1, p)) * s2 / (8.0 * pi * area);
    res[2] = (eta3 * omega - zeta * dot3(gamma2, p)) * s3 / (8.0 * pi * area);
    res[3] = (eta1 * omega - zeta * dot3(gamma3, p)) * s1 / (8.0 * pi * area);

    return res
end

#assemble the matrix D, G and K1, and K2
function assemble_matirx_DGK1K2(demag::DemagFEM, sim::MicroSimFEM)
    mesh = sim.mesh
    D = demag.D
    G = demag.G
    K1 = spzeros(mesh.number_nodes, mesh.number_nodes)
    K2 = spzeros(mesh.number_nodes, mesh.number_nodes)
    Ms = sim.Ms
    map_g2b = mesh.map_g2b

    for c = 1 : mesh.number_cells

        k1 = mesh.cell_verts[1, c]
        k2 = mesh.cell_verts[2, c]
        k3 = mesh.cell_verts[3, c]
        k4 = mesh.cell_verts[4, c]

        v = (k1, k2, k3, k4)

        bcd = compute_bcd12(mesh.coordinates[:, k1], 
                        mesh.coordinates[:, k2],
                        mesh.coordinates[:, k3],
                        mesh.coordinates[:, k4])

        for alpha = 1:4, beta = 1:4
            i = 3 * (alpha-1)
            j = 3 * (beta-1)

            fa = mesh.volumes[c]*(bcd[i+1] * bcd[j+1] + bcd[i + 2] * bcd[j + 2] + bcd[i + 3] * bcd[j + 3]);

            K1[v[alpha], v[beta]] += fa

            if map_g2b[v[alpha]] < 0
                K2[v[alpha], v[beta]] += fa
            end

            for k=1:3
                dg1 = bcd[j + k] * mesh.volumes[c] / 4.0;
                dg2 = dg1 * Ms[c];

                D[v[beta], 3*(v[alpha]-1) + k] += dg2

                G[3*(v[alpha]-1) + k, v[beta]] += dg1
            end
        end
    end

    for i = 1:mesh.number_nodes_bnd
        j = mesh.map_b2g[i]
        K2[j,j] = 1.0
    end

    demag.K1 = factorize(K1) #reuse the factorization to speed up the calculation
    demag.K2 = factorize(K2)
end

function ComputeVertBSA(mesh::FEMesh)

  vertbsa = zeros(mesh.number_nodes)
  for c = 1:mesh.number_cells, j = 0:3
    
    k1 =  mesh.cell_verts[(j+0)%4 + 1, c]
    k2 =  mesh.cell_verts[(j+1)%4 + 1, c]
    k3 =  mesh.cell_verts[(j+2)%4 + 1, c]
    k4 =  mesh.cell_verts[(j+3)%4 + 1, c]

    v1 = mesh.coordinates[:, k1]
    v2 = mesh.coordinates[:, k2]
    v3 = mesh.coordinates[:, k3]
    v4 = mesh.coordinates[:, k4]
     
    omega = solid_angle_single(v1, v2, v3, v4);

    vertbsa[k1] += omega
  end

  return vertbsa;
end


function assmeble_matrix_B(demag::DemagFEM, sim::MicroSimFEM)
    mesh = sim.mesh
    B = demag.B
    map_g2b = mesh.map_g2b

    Threads.@threads for n = 1 : mesh.number_nodes

        #if the node is not at the boundary, do nothing
        if map_g2b[n] < 0
            continue
        end

        #local v = mesh.coordinates[:, n]

        for f = 1:mesh.number_faces_bnd
            local i = mesh.face_verts[1, f]
            local j = mesh.face_verts[2, f]
            local k = mesh.face_verts[3, f]

            # skip if the node belongs to the triangle
            if n == i || n == j || n == k
                continue
            end

            local be = boundary_element(mesh.coordinates[:, n], 
                                        mesh.coordinates[:, i],
                                        mesh.coordinates[:, j],
                                        mesh.coordinates[:, k])


            B[map_g2b[n], map_g2b[i]] += be[1]
            B[map_g2b[n], map_g2b[j]] += be[2]
            B[map_g2b[n], map_g2b[k]] += be[3]
        end
    end

    vertbsa = ComputeVertBSA(mesh)

    for n = 1 : mesh.number_nodes
        if map_g2b[n] > 0
            tmp = vertbsa[n] / (4.0 * pi) - 1.0
            B[map_g2b[n], map_g2b[n]] += tmp
        end
    end

end


function init_demag(sim::MicroSimFEM)

    mesh = sim.mesh
    
    demag = DemagFEM()
    
    #demag.K1 = spzeros(mesh.number_nodes, mesh.number_nodes)
    #demag.K2 = spzeros(mesh.number_nodes, mesh.number_nodes)
    demag.D = spzeros(mesh.number_nodes, 3*mesh.number_nodes)
    demag.G = spzeros(3*mesh.number_nodes, mesh.number_nodes)
  
    demag.B = zeros(mesh.number_nodes_bnd, mesh.number_nodes_bnd)
    
    demag.g1 = zeros(mesh.number_nodes)
    demag.g2 = zeros(mesh.number_nodes)
    demag.phi1 = zeros(mesh.number_nodes)
    demag.phi2 = zeros(mesh.number_nodes)
  
    demag.u1_bnd = zeros(mesh.number_nodes_bnd)
    demag.u2_bnd = zeros(mesh.number_nodes_bnd)
  
    demag.field = zeros(3*mesh.number_nodes)
    demag.energy = zeros(mesh.number_nodes)
  
    demag.matrices_assembled = false
  
    return demag
  
  end