using JuMag
using Test

function compute_field_macro_uniform(mesh, Nx, Ny, Nz; gpu=false)
    Ms = 8.6e5
    sim = Sim(mesh)
    sim.Ms[:] .= Ms

    function rand_m(i, j, k, dx, dy, dz)
      return (1, sin(i*0.2), cos(j*0.3))
    end

    init_m0(sim, (1,1,1), norm=false)
    add_demag(sim, Nx=Nx, Ny=Ny, Nz=Nz)
    if gpu
        JuMag.compute_fields_to_gpu(sim, sim.spin, 0.0)
    else
        JuMag.effective_field(sim, sim.spin, 0.0)
    end
    return sim.field/Ms
end


function compute_field_macro(mesh, Nx, Ny, Nz; gpu=false)
    Ms = 8.6e5
    sim = Sim(mesh)
    sim.Ms[:] .= Ms

    function rand_m(i, j, k, dx, dy, dz)
        I = i%3
        J = j%3
        return (1, sin(I*0.8), cos(J*0.7))
    end

    init_m0(sim, rand_m)
    add_demag(sim, Nx=Nx, Ny=Ny, Nz=Nz)
    if gpu
        JuMag.compute_fields_to_gpu(sim, sim.spin, 0.0)
    else
        JuMag.effective_field(sim, sim.spin, 0.0)
    end
    return sim.field/Ms
end

function test_field_macro(;gpu=false)
    if gpu
        mesh1 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=1)
        mesh2 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=3, nz=1)
        mesh3 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=101, ny=101, nz=1)
        mesh4 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=33, ny=33, nz=1)
    else
        mesh1 =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=1)
        mesh2 =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=3, nz=1)
        mesh3 =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=101, ny=101, nz=1)
        mesh4 =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=33, ny=33, nz=1)
    end
    #=
    f1 = compute_field_macro_uniform(mesh1, 500, 0, 0, gpu=gpu)  #1d case
    println(f1)
    @test isapprox(f1[1], 0, atol=1e-6)
    @test isapprox(f1[2], -0.5, atol=1e-6)
    @test isapprox(f1[3], -0.5, atol=1e-6)

    f2 = compute_field_macro_uniform(mesh3, 5, 5, 0, gpu=gpu)  #2d case
    field = reshape(f2, 3, 101, 101, 1)
    f2b = field[:, 51, 51, 1]
    println(f2b)
    @test isapprox(f2b[1], 0, atol=1e-3)
    @test isapprox(f2b[2], 0, atol=1e-3)
    @test isapprox(f2b[3], -1, atol=1e-3)

    f3 = compute_field_macro_uniform(mesh1, 31, 31, 0, gpu=gpu) #self-consistency
    f4 = compute_field_macro_uniform(mesh2, 10, 10, 0, gpu=gpu)
    field = reshape(f4, 3, 3, 3, 1)
    f3b = field[:, 2, 2, 1]
    println(f3, f3b)
    @test isapprox(f3, f3b, atol=1e-10)
    =#
    f5 = compute_field_macro_uniform(mesh2, 5, 5, 0, gpu=gpu) #self-consistency
    f6 = compute_field_macro_uniform(mesh4, 0, 0, 0, gpu=gpu)
    field = reshape(f6, 3, 33, 33, 1)
    f6b = field[1:3, 16:18, 16:18, 1]
    field = reshape(f6b, 27)
    println(size(f5), size(f5), typeof(f5), typeof(f6b))
    println(field1.-f6b)
    #diff = field1[1:3, 1:3, 1:3, 1] - f6b
    #println(diff)
    @info(f6b)

end

@testset "Macro Boundary" begin
    test_field_macro()
end

if JuMag._cuda_available.x
  JuMag.cuda_using_double()
  @testset "Macro Boundary GPU" begin
      test_field_macro(gpu=true)
  end
end
