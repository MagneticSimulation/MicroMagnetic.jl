@inline function distance(x1::Array{Float64,1}, x2::Array{Float64,1})
    x = x1 .- x2
    return norm2(x)
end

@inline function kronecker_delta(x::Int64, y::Int64)
    if x == y
        return 1
    end
    return 0
end

#
# helper function, compute the det of the given three vectors
# for equation (9), equals 6*a1*Ve
#
function det3(x::Array{Float64,1}, y::Array{Float64,1}, z::Array{Float64,1})
    d = x[1] * y[2] * z[3] + x[2] * y[3] * z[1] + x[3] * y[1] * z[2] - x[1] * y[3] * z[2] -
        x[2] * y[1] * z[3] - x[3] * y[2] * z[1]
    return d
end

function norm2(X::Array{Float64,1})
    v = 0.0
    for x in X
        v += x * x
    end
    return sqrt(v)
end

function dot3(x::Array{Float64,1}, y::Array{Float64,1})
    return x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
end

function triangle_area(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1})
    a = x2 .- x1
    b = x3 .- x2
    c = x1 .- x3
    s1 = norm2(a)
    s2 = norm2(b)
    s3 = norm2(c)
    s = (s1 + s2 + s3) / 2.0
    area = sqrt(s * (s - s1) * (s - s2) * (s - s3))
    return area
end

function cross3(x::Array{Float64,1}, y::Array{Float64,1})
    return [-x[3] * y[2] + x[2] * y[3], x[3] * y[1] - x[1] * y[3],
            -x[2] * y[1] + x[1] * y[2]]
end



#compute the volume of a tetrahedron
function tetrahedron_volume(x1::Array{Float64,1}, x2::Array{Float64,1},
                            x3::Array{Float64,1}, x4::Array{Float64,1})
    v = det3(x2, x3, x4) - det3(x1, x3, x4) + det3(x1, x2, x4) - det3(x1, x2, x3)
    return v / 6.0
end

function tet_volume(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1},
                    x4::Array{Float64,1})
    v = det3(x2 .- x1, x3 .- x1, x4 .- x1)
    return v / 6.0
end

#
# equation (10), (11) and (12), note the factor of 6*Ve
# compute the coefficients b c and d of a tetrahedron
# return (b1,c1,d1)
#
function bcd_negative(x2::Array{Float64,1}, x3::Array{Float64,1}, x4::Array{Float64,1})
    b1 = x2[3] * x3[2] - x2[2] * x3[3] - x2[3] * x4[2] + x3[3] * x4[2] + x2[2] * x4[3] -
         x3[2] * x4[3]

    c1 = -x2[3] * x3[1] + x2[1] * x3[3] + x2[3] * x4[1] - x3[3] * x4[1] - x2[1] * x4[3] +
         x3[1] * x4[3]

    d1 = x2[2] * x3[1] - x2[1] * x3[2] - x2[2] * x4[1] + x3[2] * x4[1] + x2[1] * x4[2] -
         x3[1] * x4[2]

    return (b1, c1, d1)
end

#
# compute alpha \times k \times beta, where alpha,beta=1,2,3,4; k=1,2,3;
# return b1 c1 d1 b2 c2 d2 b3 c3 d3 b4 c4 d4
# 
function compute_bcd12(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1},
                       x4::Array{Float64,1})
    J = 1.0 / (6 * tet_volume(x1, x2, x3, x4))

    (b1, c1, d1) = bcd_negative(x2, x3, x4) .* J
    (b2, c2, d2) = bcd_negative(x1, x3, x4) .* (-J)
    (b3, c3, d3) = bcd_negative(x2, x1, x4) .* (-J)
    (b4, c4, d4) = bcd_negative(x2, x3, x1) .* (-J)

    return (b1, c1, d1, b2, c2, d2, b3, c3, d3, b4, c4, d4)
end

# compute the coefficients b1 c1 and d1 of a tetrahedron
function compute_bcd1(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1},
                      x4::Array{Float64,1})
    J = 1.0 / (6 * tet_volume(x1, x2, x3, x4))
    return bcd_negative(x2, x3, x4) .* J
end

function compute_correction_triangle(dx::Vector{Float64}, dy::Vector{Float64},
                                     dz::Vector{Float64}, sa::Float64, sb::Float64,
                                     sc::Float64)::Float64
    r1 = distance(dx, dy)
    r2 = distance(dx, dz)
    r3 = distance(dy, dz)

    fa = det2(dx, dy, dz)

    fb = (sb - sc) * (r1 - r2) / (2 * r3 * r3)
    fc = (sb + sc + 2 * sa) / (4.0 * r3) + (r2^2 - r1^2) * (sb - sc) / (4.0 * r3^3)
    fd = log(r1 + r2 + r3) - log(r1 + r2 - r3)

    return fa * (fb + fc * fd)
end

function compute_correction_triangle_A(dx::Vector{Float64}, dy::Vector{Float64},
                                       dz::Vector{Float64}, sa::Float64)::Float64
    r1 = distance(dx, dy)
    r2 = distance(dx, dz)
    r3 = distance(dy, dz)

    fa = det2(dx, dy, dz)

    fc = (sa) / (2.0 * r3)
    fd = log(r1 + r2 + r3) - log(r1 + r2 - r3)

    return fa * fc * fd
end

# compute the mesh volumes of given mesh
function compute_mesh_volume!(mesh::FEMesh)
    volumes = zeros(mesh.number_cells)

    for i in 1:(mesh.number_cells)
        k1 = mesh.cell_verts[1, i]
        k2 = mesh.cell_verts[2, i]
        k3 = mesh.cell_verts[3, i]
        k4 = mesh.cell_verts[4, i]

        v = tetrahedron_volume(mesh.coordinates[:, k1], mesh.coordinates[:, k2],
                               mesh.coordinates[:, k3], mesh.coordinates[:, k4])

        volumes[i] = abs(v)
    end

    return mesh.volumes = volumes
end

function compute_L_inv_neg!(mesh::FEMesh)
    nodal_volumes = zeros(mesh.number_nodes)
    cv = mesh.cell_verts
    for i in 1:(mesh.number_cells), k in 1:4
        nodal_volumes[cv[k, i]] += mesh.volumes[i] / 4.0
    end

    L_inv_neg = zeros(3 * mesh.number_nodes)

    for i in 0:(mesh.number_nodes - 1)
        L_inv_neg[3 * i + 1] = -1.0 / nodal_volumes[i + 1]
        L_inv_neg[3 * i + 2] = -1.0 / nodal_volumes[i + 1]
        L_inv_neg[3 * i + 3] = -1.0 / nodal_volumes[i + 1]
    end

    return mesh.L_inv_neg = kernel_array(L_inv_neg)
end

function compute_L_Ms!(L_mu::AbstractArray{T,1}, mesh::FEMesh, Ms::Array{T,1}) where {T<:AbstractFloat}
    nodal_Ms = zeros(mesh.number_nodes)
    cv = mesh.cell_verts

    for k in 1:4, i in 1:(mesh.number_cells)
        nodal_Ms[cv[k, i]] += mesh.volumes[i] * Ms[i] / 4.0
    end

    L_mu_tmp = Array(L_mu)
    for i in 0:(mesh.number_nodes - 1)
        L_mu_tmp[3 * i + 1] = nodal_Ms[i + 1]
        L_mu_tmp[3 * i + 2] = nodal_Ms[i + 1]
        L_mu_tmp[3 * i + 3] = nodal_Ms[i + 1]
    end
    copyto!(L_mu, L_mu_tmp)
end

function compute_normal(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1})
    v1 = x2 .- x1
    v2 = x3 .- x1
    v = cross3(v1, v2)
    return v / norm2(v)
end


function init_vector!(v::Array{T,1}, mesh::FEMesh, init::Function) where {T<:AbstractFloat}
    N = mesh.number_nodes
    b = reshape(v, 3, N)

    for i in 1:N
        x, y, z = mesh.coordinates[:, i]
        b[:, i] .= init(x, y, z)
    end

    if NaN in v
        error("NaN is given by the input function.")
    end
    return nothing
end

function init_scalar_nodes!(v::Array{T,1}, mesh::FEMesh, init::Function) where {T<:AbstractFloat}
    N = mesh.number_nodes

    for i in 1:N
        x, y, z = mesh.coordinates[:, i]
        v[i] = init(x, y, z)
    end

    if NaN in v
        error("NaN is given by the input function.")
    end
    return nothing
end


function init_vector!(v::AbstractArray{T,1}, mesh::FEMesh, init_fun) where {T<:AbstractFloat}
    init_v = zeros(T, 3 * mesh.number_nodes)
    init_vector!(init_v, mesh, init_fun)
    copyto!(v, init_v)
    return true
end

function init_vector!(v::Array{T,1}, mesh::FEMesh,
                      init::Tuple{Real,Real,Real}) where {T<:AbstractFloat}
    #N = length(v)
    N = mesh.number_nodes
    b = reshape(v, 3, N)
    b[1, :] .= init[1]
    b[2, :] .= init[2]
    b[3, :] .= init[3]
    return nothing
end

function init_vector!(v::Array{T,1}, mesh::FEMesh,
                      init::Array{T,1}) where {T<:AbstractFloat}
    v[:] = init[:]
    return nothing
end

function init_scalar!(v::Array{T,1}, mesh::FEMesh, init::Number) where {T<:AbstractFloat}
    for i in 1:(mesh.number_cells)
        v[i] = init
    end
    return nothing
end

function init_scalar!(v::Array{T,1}, mesh::FEMesh, init::Function) where {T<:AbstractFloat}
    for i in 1:(mesh.number_cells)
        material_id = mesh.material_ids[i]
        v[i] = init(i, material_id)
    end
    return nothing
end

function init_scalar!(v::Array{T,1}, mesh::FEMesh,
                      init::Array{T,1}) where {T<:AbstractFloat}
    v[:] = init[:]
    return nothing
end

function assemble_mass_matirx(mesh)
    N = mesh.number_nodes
    M = spzeros((3 * N, 3 * N))

    for c in 1:(mesh.number_cells)
        k1 = mesh.cell_verts[1, c]
        k2 = mesh.cell_verts[2, c]
        k3 = mesh.cell_verts[3, c]
        k4 = mesh.cell_verts[4, c]

        v = ((k1, k2, k3, k4) .- 1) .* 3

        for alpha in 1:4, beta in 1:4
            coef = mesh.volumes[c] * (1 + kronecker_delta(alpha, beta)) / (20.0)

            for i in 1:3
                M[v[alpha] + i, v[beta] + i] += coef
            end
        end
    end

    return M
end

function compute_character_length(mesh)
    v = sum(mesh.volumes) / length(mesh.volumes)

    return v^(1 / 3)
end
