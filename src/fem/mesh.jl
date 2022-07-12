mutable struct FEMesh <: Mesh
    number_nodes::Int64  #total node number
    number_cells::Int64  #total cell number
    number_faces_bnd::Int64 #total surface number
    coordinates::Array{Float64, 2}  #coordinates array
    cell_verts::Array{Int64, 2} 
    face_verts::Array{Int64, 2} 
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

    for i = 1:N
        x = split(readline(io))
        #print(x)
        @assert i == parse(Int64, x[1])
        mesh.cell_verts[1, i] = parse(Int64, x[2])
        mesh.cell_verts[2, i] = parse(Int64, x[3])
        mesh.cell_verts[3, i] = parse(Int64, x[4])
        mesh.cell_verts[4, i] = parse(Int64, x[5])
    end

    N = parse(Int64, readline(io)) #finally it sould be the total surface number
    mesh.number_faces_bnd = N
    mesh.face_verts = zeros(Int64, 3, N)

    for i = 1:N
        x = split(readline(io))
        j = parse(Int64, x[1])
        mesh.face_verts[1, j] = parse(Int64, x[2])
        mesh.face_verts[2, j] = parse(Int64, x[3])
        mesh.face_verts[3, j] = parse(Int64, x[4])
    end
    return mesh
end

function FEMesh(fname::String)
    return load_mesh_netgen_neutral(fname)
end