export FEMesh

mutable struct FEMesh <: Mesh
    number_nodes::Int64  #total node number
    number_cells::Int64  #total cell number
    number_faces_bnd::Int64 #total surface number
    number_nodes_bnd::Int64 #total node number at surface
    unit_length::Float64
    coordinates::Array{Float64,2}  #coordinates array
    cell_verts::Array{Int32,2}
    face_verts::Array{Int32,2}
    map_g2b::AbstractArray{Int32,1}
    map_b2g::AbstractArray{Int32,1}
    volumes::Array{Float64,1}
    L_inv_neg::AbstractArray{Float64,1}
    region_ids::Array{Int64,1}
    surface_ids::Array{Int64,1}
    surface_normals::Array{Float64,2}
    FEMesh() = new()
end

function UnitTetrahedronMesh()
    mesh = FEMesh()
    mesh.number_nodes = 4
    mesh.coordinates = [0.0 1 0 0;
                        0.0 0 1 0;
                        0.0 0 0 1]

    mesh.number_cells = 1
    mesh.cell_verts = reshape([1; 2; 3; 4], 4, 1)

    mesh.number_faces_bnd = 4
    mesh.face_verts = [1 1 1 2;
                       2 2 3 3;
                       3 4 4 4]

    mesh.unit_length = 1.0
    return mesh
end

function load_mesh_netgen_neutral(fname::String)
    io = open(fname, "r")

    mesh = FEMesh()

    N = parse(Int64, readline(io)) #the first line should be the total number of nodes
    mesh.number_nodes = N
    mesh.coordinates = zeros(Float64, 3, N)
    for i in 1:N
        x = split(readline(io))
        #print(x)
        mesh.coordinates[1, i] = parse(Float64, x[1])
        mesh.coordinates[2, i] = parse(Float64, x[2])
        mesh.coordinates[3, i] = parse(Float64, x[3])
    end

    N = parse(Int64, readline(io)) #then it sould be the total number of cells
    mesh.number_cells = N
    mesh.cell_verts = zeros(Float64, 4, N)
    mesh.region_ids = zeros(Integer, N)

    for i in 1:N
        x = split(readline(io))
        mesh.region_ids[i] = parse(Int64, x[1])
        mesh.cell_verts[1, i] = parse(Int64, x[2])
        mesh.cell_verts[2, i] = parse(Int64, x[3])
        mesh.cell_verts[3, i] = parse(Int64, x[4])
        mesh.cell_verts[4, i] = parse(Int64, x[5])
    end

    N = parse(Int64, readline(io)) #finally it sould be the total surface number
    mesh.number_faces_bnd = N
    mesh.face_verts = zeros(Int64, 3, N)
    mesh.surface_ids = zeros(Integer, N)

    for i in 1:N
        x = split(readline(io))
        mesh.surface_ids[i] = parse(Int64, x[1])
        mesh.face_verts[1, i] = parse(Int64, x[2])
        mesh.face_verts[2, i] = parse(Int64, x[3])
        mesh.face_verts[3, i] = parse(Int64, x[4])
    end
    return mesh
end

function build_boundary_maps!(mesh::FEMesh)
    map_g2b = zeros(Int32, mesh.number_nodes)
    map_g2b .= -1

    at_bounary = zeros(Bool, mesh.number_nodes)
    for i in mesh.face_verts
        at_bounary[i] = true
    end

    number_nodes_bnd = 0
    for i in 1:(mesh.number_nodes)
        if at_bounary[i]
            number_nodes_bnd += 1
            map_g2b[i] = number_nodes_bnd
        end
    end
    mesh.number_nodes_bnd = number_nodes_bnd
    mesh.map_g2b = kernel_array(map_g2b)

    map_b2g = zeros(Int32, mesh.number_nodes_bnd)
    for i in 1:(mesh.number_nodes)
        if map_g2b[i] > 0
            map_b2g[map_g2b[i]] = i
        end
    end
    return mesh.map_b2g = kernel_array(map_b2g)
end

function compute_surface_normals!(mesh::FEMesh)
    surface_normals = zeros(Float64, 3, mesh.number_faces_bnd)
    for i in 1:(mesh.number_faces_bnd)
        v1 = mesh.coordinates[:, mesh.face_verts[1, i]]
        v2 = mesh.coordinates[:, mesh.face_verts[2, i]]
        v3 = mesh.coordinates[:, mesh.face_verts[3, i]]
        n = cross3(v2 - v1, v3 - v1)
        n = n / norm2(n)
        surface_normals[:, i] = n
    end
    mesh.surface_normals = surface_normals
    return mesh.surface_normals
end

"""
Create a FEMesh from the given netgen neutral file.
"""
function FEMesh(fname::String; unit_length=1e-9)
    mesh = load_mesh_netgen_neutral(fname)
    mesh.unit_length = unit_length
    compute_mesh_volume!(mesh)
    compute_L_inv_neg!(mesh)
    build_boundary_maps!(mesh)
    return mesh
end

function export_mesh(mesh, filename)
    f = open(filename, "w")

    # write the number of nodes
    write(f, string(mesh.number_nodes, "\n"))
    for i in 1:(mesh.number_nodes)
        x = mesh.coordinates[:, i]
        write(f, string(x[1], " ", x[2], " ", x[3], "\n"))
    end

    # write the number of cells (tetrahedra here)
    write(f, string(mesh.number_cells, "\n"))
    for i in 1:(mesh.number_cells)
        x = mesh.cell_verts[:, i]
        material = mesh.region_ids[i]
        write(f, string(material, " ", x[1], " ", x[2], " ", x[3], " ", x[4], "\n"))
    end

    write(f, string(mesh.number_faces_bnd, "\n"))
    for i in 1:(mesh.number_faces_bnd)
        x = mesh.face_verts[:, i]
        write(f, string(mesh.surface_ids[i], " ", x[1], " ", x[2], " ", x[3], "\n"))
    end

    return close(f)
end

function merge_mesh(base_mesh, mesh)
    base_mesh.number_cells += mesh.number_cells
    cell_verts = mesh.cell_verts .+ base_mesh.number_nodes
    base_mesh.cell_verts = hcat(base_mesh.cell_verts, cell_verts)

    number_material = length(Set(base_mesh.region_ids)) + length(Set(mesh.region_ids))
    #base_mesh.number_material += mesh.number_material  
    region_ids = mesh.region_ids .+ (number_material - 1)
    base_mesh.region_ids = vcat(base_mesh.region_ids, region_ids)

    base_mesh.number_faces_bnd += mesh.number_faces_bnd
    face_verts = mesh.face_verts .+ base_mesh.number_nodes
    base_mesh.face_verts = hcat(base_mesh.face_verts, face_verts)
    surface_ids = mesh.surface_ids .+ length(Set(base_mesh.surface_ids))
    base_mesh.surface_ids = vcat(base_mesh.surface_ids, surface_ids)

    base_mesh.number_nodes += mesh.number_nodes
    base_mesh.coordinates = hcat(base_mesh.coordinates, mesh.coordinates)
    return base_mesh
end
