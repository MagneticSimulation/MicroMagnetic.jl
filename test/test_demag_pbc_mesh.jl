#Tests for the mesh-driven demag PBC dispatch: add_demag(sim) without PBC
#kwargs follows the mesh periodicity (FDMesh(...; pbc=...)) and must agree
#with the explicit solver choice; conflicting explicit kwargs warn.

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
    #1d: mesh pbc="y" dispatches to pbc1d on y
    nx, ny, nz = 6, 5, 3
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz, pbc="y")
    m = mesh_random_m(nx, ny, nz; seed=1)
    h1 = field_of(mesh, m)
    h2 = field_of(mesh, m; pbc1d=true, Ny=2)
    scale = maximum(abs, h2)
    tol = eltype(h2) == Float32 ? 1e-4 : 1e-9
    @test maximum(abs, h1 - h2) <= tol * scale

    #2d: mesh pbc="xz" dispatches to pbc2d on x,z
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz, pbc="xz")
    m = mesh_random_m(nx, ny, nz; seed=2)
    h1 = field_of(mesh, m)
    h2 = field_of(mesh, m; pbc2d=true, Nx=2, Nz=2)
    scale = maximum(abs, h2)
    @test maximum(abs, h1 - h2) <= tol * scale

    #3d: mesh pbc="xyz" dispatches to pbc3d
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=4, ny=4, nz=4, pbc="xyz")
    m = mesh_random_m(4, 4, 4; seed=3)
    h1 = field_of(mesh, m)
    h2 = field_of(mesh, m; pbc3d=true)
    scale = maximum(abs, h2)
    @test maximum(abs, h1 - h2) <= tol * scale

    #open mesh: the bare call gives the open demag
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz)
    m = mesh_random_m(nx, ny, nz; seed=4)
    h1 = field_of(mesh, m)
    h2 = field_of(mesh, m; Nx=0, Ny=0, Nz=0)
    @test maximum(abs, h1 - h2) <= tol * maximum(abs, h2)

    #conflicting explicit kwargs warn (mesh xy, demag 1d on x)
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz, pbc="xy")
    m = mesh_random_m(nx, ny, nz; seed=5)
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, m; norm=false)
    @test_logs (:warn,) match_mode=:any add_demag(sim; pbc1d=true, Nx=2)
end

@using_gpu()
test_functions("DemagPBCMesh", test_demag_pbc_mesh)
