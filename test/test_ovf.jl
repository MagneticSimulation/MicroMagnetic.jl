using JuMag
using Test

mesh = FDMesh(nx = 2, ny = 1, nz = 1, dx = 1e-9, dy = 1e-9, dz = 1e-9)
sim = Sim(mesh, name = "test_ovf1")
set_Ms(sim, 8.0e5)
init_m0(sim, (0.6,0.8,0))
save_ovf(sim, "test_ovf_bin", dataformat = "binary")
save_ovf(sim, "test_ovf_text", dataformat = "text")

function test_read_ovf(ovf_name)
    sim = Sim(mesh, name = "test_ovf3")
    set_Ms(sim, 8.0e5)
    read_ovf(sim, ovf_name)
    println(sim.spin)
    @test sim.spin[1] == 0.6
    @test sim.spin[2] == 0.8
    @test sim.spin[3] == 0.0
    @test sim.prespin[1] == 0.6
    @test sim.prespin[2] == 0.8
    @test sim.prespin[3] == 0.0
end

test_read_ovf("test_ovf_bin")
test_read_ovf("test_ovf_text")
#JuMag.cuda_using_double(true)
