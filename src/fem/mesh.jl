mutable struct FEMesh <: Mesh
    number_nodes::Int64  #total node number
    number_cells::Int64  #total cell number
    number_faces_bnd::Int64 #total surface number
    number_nodes_bnd::Int64 #total node number at surface
    unit_length::Float64
    coordinates::Array{Float64, 2}  #coordinates array
    cell_verts::Array{Int64, 2} 
    face_verts::Array{Int64, 2}
    map_g2b::Array{Int64, 1}
    map_b2g::Array{Int64, 1}
    volumes::Array{Float64, 1}
    L_inv_neg::Array{Float64, 1}
    material_ids::Array{Int64, 1}
    surface_ids::Array{Int64, 1}
    FEMesh() = new()
end

function UnitTetrahedronMesh()
    mesh = FEMesh()
    mesh.number_nodes = 4
    mesh.coordinates = [0.0 1 0 0;
                   0.0 0 1 0;
                   0.0 0 0 1];
    
    mesh.number_cells = 1
    mesh.cell_verts = reshape([1;2;3;4], 4, 1)
    
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
    for i = 1:N
        x = split(readline(io))
        #print(x)
        mesh.coordinates[1,i] = parse(Float64, x[1])
        mesh.coordinates[2,i] = parse(Float64, x[2])
        mesh.coordinates[3,i] = parse(Float64, x[3])
    end

    N = parse(Int64, readline(io)) #then it sould be the total number of cells
    mesh.number_cells = N
    mesh.cell_verts = zeros(Float64, 4, N)
    mesh.material_ids = zeros(Integer, N)

    for i = 1:N
        x = split(readline(io))
        mesh.material_ids[i] = parse(Int64, x[1])
        mesh.cell_verts[1, i] = parse(Int64, x[2])
        mesh.cell_verts[2, i] = parse(Int64, x[3])
        mesh.cell_verts[3, i] = parse(Int64, x[4])
        mesh.cell_verts[4, i] = parse(Int64, x[5])
    end

    N = parse(Int64, readline(io)) #finally it sould be the total surface number
    mesh.number_faces_bnd = N
    mesh.face_verts = zeros(Int64, 3, N)
    mesh.surface_ids = zeros(Integer, N)

    for i = 1:N
        x = split(readline(io))
        mesh.surface_ids[i] = parse(Int64, x[1])
        mesh.face_verts[1, i] = parse(Int64, x[2])
        mesh.face_verts[2, i] = parse(Int64, x[3])
        mesh.face_verts[3, i] = parse(Int64, x[4])
    end
    return mesh
end


function build_boundary_maps!(mesh::FEMesh)
    map_g2b = zeros(Int64, mesh.number_nodes)
    map_g2b .= -1

    at_bounary = zeros(Bool, mesh.number_nodes)
    for i in mesh.face_verts
        at_bounary[i] = true
    end

    number_nodes_bnd = 0
    for i = 1:mesh.number_nodes
        if at_bounary[i]
            number_nodes_bnd += 1
            map_g2b[i] = number_nodes_bnd
        end
    end
    mesh.number_nodes_bnd = number_nodes_bnd
    mesh.map_g2b = map_g2b

    map_b2g = zeros(Int64, mesh.number_nodes_bnd)
    for i = 1:mesh.number_nodes
        if map_g2b[i]>0
            map_b2g[map_g2b[i]] = i
        end
    end
    mesh.map_b2g = map_g2b

end

function FEMesh(fname::String; unit_length=1e-9)
    mesh = load_mesh_netgen_neutral(fname)
    mesh.unit_length = unit_length
    compute_mesh_volume!(mesh)
    compute_L_inv_neg!(mesh)
    build_boundary_maps!(mesh)
    return mesh
end


