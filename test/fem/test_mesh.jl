using JuMag

function test_tetrahedron_mesh()
    mesh1 = UnitTetrahedronMesh()
    filepath = joinpath(@__DIR__, "tet.neutral")
    #print(filepath)
    mesh2 = FEMesh(filepath)
    @assert maximum(mesh1.coordinates-mesh2.coordinates)==0.0
end

function test_read_neutral()
    filepath = joinpath(@__DIR__, "octa.neutral")
    mesh = FEMesh(filepath)
    @assert mesh.number_nodes == 7
    @assert mesh.number_cells == 8
    @assert mesh.number_faces_bnd == 8
    #@assert mesh.volumes[4] == 1.0/6
end

test_tetrahedron_mesh()

test_read_neutral()

