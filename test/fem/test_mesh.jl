using JuMag
using LinearAlgebra

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

function test_extract_normal_axes_by_maximum_area()
    filepath = joinpath(@__DIR__, "two_hex.mesh")
    mesh = FEMesh(filepath)
    normals = JuMag.extract_normal_axes_by_maximum_area(mesh)
    expected_n1 = [-0.8707639652559477, -0.4888975853320629, 0.05243346134117942]
    expected_n2 = [0.0022083576995702444, -0.02397243455089625, 0.9997101807713943]
    print(dot(cross(normals[:, 1], expected_n1),expected_n1))
    @assert dot(cross(normals[:, 1], expected_n1),expected_n1) < 1e-15
    @assert dot(cross(normals[:, 2], expected_n2),expected_n2) < 1e-15
end

test_tetrahedron_mesh()

test_read_neutral()

test_extract_normal_axes_by_maximum_area()

