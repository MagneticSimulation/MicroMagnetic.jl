#Tests for the mesh-driven demag PBC dispatch: the bare add_demag follows the
#mesh periodicity; macroPBC=true overrides it with the legacy truncated-image
#method; a macro periodicity that differs from the mesh warns.

using MicroMagnetic
using Random
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function mesh_random_m(nx, ny, nz; seed)
    Random.seed!(seed)
    m = randn(3 * nx * ny * nz)
    for c in 1:3
        idx = c:3:length(m)
        m[idx] .-= sum(m[idx]) / length(idx)
    end
    return m
end

function field_of(mesh, m; kwargs...)
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, m; norm=false)
    add_demag(sim; kwargs...)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    return Array(sim.field)
end

function test_demag_pbc_mesh()
    #the open mesh: the bare call and macro=true (no counts) both give the open demag
    nx, ny, nz = 6, 5, 3
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz)
    m = mesh_random_m(nx, ny, nz; seed=4)
    h1 = field_of(mesh, m)
    h2 = field_of(mesh, m; macroPBC=true)
    @test maximum(abs, h1 - h2) <= 1e-9 * maximum(abs, h2)

    #macroPBC=true overrides the mesh pbc: the mesh pbc="y" + macroPBC, Ny=2
    #gives the macro field -- close to, but distinct from, the true-pbc1d field
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz, pbc="y")
    m = mesh_random_m(nx, ny, nz; seed=1)
    h_true = field_of(mesh, m)                       #pbc1d, Ic=4
    h_mac = field_of(mesh, m; macroPBC=true, Ny=2)
    scale = maximum(abs, h_true)
    d = maximum(abs, h_mac - h_true)
    Tf = eltype(h_true)
    @test d > 0                                      #the override took effect
    @test d <= (Tf == Float32 ? 1e-1 : 1e-2) * scale #... and stays near the true field

    #a macro periodicity that differs from the mesh warns
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz, pbc="xy")
    m = mesh_random_m(nx, ny, nz; seed=5)
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, m; norm=false)
    @test_logs (:warn,) match_mode=:any add_demag(sim; macroPBC=true, Nx=2)
end

@using_gpu()
test_functions("DemagPBCMesh", test_demag_pbc_mesh)
