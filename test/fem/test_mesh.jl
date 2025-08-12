using MicroMagnetic
using LinearAlgebra

function test_read_neutral()
    filepath = joinpath(@__DIR__, "meshes/octa.mesh")
    mesh = FEMesh(filepath)
    @assert mesh.number_nodes == 7
    @assert mesh.number_cells == 8
    @assert mesh.number_faces_bnd == 8
    #@assert mesh.volumes[4] == 1.0/6
end

test_read_neutral()