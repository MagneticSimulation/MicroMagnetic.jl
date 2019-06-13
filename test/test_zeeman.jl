using JuMag
using Test

function fun_x(t)
    return sin(1e9*t)
end

function fun_y(t)
    return cos(1e9*t)
end

function fun_z(t)
    return 1.0
end

function test_zeeman(;gpu=false)
    if gpu
        mesh =  FDMeshGPU(dx=2e-9, nx=3, ny=2, nz=1, pbc="x")
    else
        mesh =  FDMesh(dx=2e-9, nx=3, ny=2, nz=1, pbc="x")
    end

    sim = Sim(mesh, name="test_zeeman")
    sim.driver.alpha = 0.01
    sim.driver.gamma = 2.21e5

    init_m0(sim, (1,1,1), norm=false)

    set_Ms(sim, 8.6e5)

    z1 = add_zeeman(sim, (1,2,2e3))
    z2 = add_zeeman(sim, (1e3,1e4,1e5), (fun_x, fun_y, fun_z))

    if gpu
        JuMag.compute_fields_to_gpu(sim, sim.spin, 1.23e-11)
    else
       JuMag.effective_field(sim, sim.spin, 1.23e-11)
    end

    @test z1.field[1] == 1.0
    @test z1.field[2] == 2.0
    @test z1.field[3] == 2e3

    @test z2.field[1] == 1e3*sin(1.23e-2)
    @test z2.field[2] == 1e4*cos(1.23e-2)
    @test z2.field[3] == 1e5
    print(z2.field)
    return nothing
end

test_zeeman()
