using JuMag
using Test

function compute_field_macro(mesh, Nx, Ny, Nz)
    Ms = 8.6e5
    sim = Sim(mesh)
    sim.Ms[:] .= Ms

    init_m0(sim, (1,1,1), norm=false)
    add_demag(sim, Nx=Nx, Ny=Ny, Nz=Nz)

    JuMag.effective_field(sim, sim.spin, 0.0)

    return sim.field/Ms
end

function test_field_macro(;using_gpu=false)
    if using_gpu
        mesh1 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=1)
        mesh2 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=3, nz=1)
    else
        mesh1 =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=1)
        mesh2 =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=3, nz=1)
    end

    f1 = compute_field_macro(mesh1, 500, 0, 0)
    println(f1)
    @test isapprox(f1[1], 0, atol=1e-6)
    @test isapprox(f1[2], -0.5, atol=1e-6)
    @test isapprox(f1[3], -0.5, atol=1e-6)

    f2 = compute_field_macro(mesh1, 500, 500, 0)
    println(f2)
    @test isapprox(f2[1], 0, atol=1e-3)
    @test isapprox(f2[2], 0, atol=1e-3)
    @test isapprox(f2[3], -1, atol=1e-3)

    f3 = compute_field_macro(mesh1, 31, 31, 0)
    f4 = compute_field_macro(mesh2, 10, 10, 0)
    field = reshape(f4, 3, 3, 3, 1)
    f3b = field[:, 2, 2, 1]
    println(f3, f3b)
    @test isapprox(f3, f3b, atol=1e-11)

end

test_field_macro()
