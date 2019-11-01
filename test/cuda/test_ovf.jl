using JuMag
using Test
JuMag.cuda_using_double(true)
mesh = FDMeshGPU(nx = 2, ny = 1, nz = 1, dx = 1e-9, dy = 1e-9, dz = 1e-9)
sim = Sim(mesh, name = "test_ovf1")
set_Ms(sim, 8.0e5)
init_m0(sim, (0.6,0.8,0))
save_ovf(sim, "test_ovf", dataformat = "binary")

expected = [0.6, 0.8, 0, 0.6, 0.8, 0.0]

JuMag.cuda_using_double(false)
sim = Sim(mesh, name = "test_ovf2")
set_Ms(sim, 8.0e5)
init_m0(sim, (0.6,0.8,0))
save_ovf(sim, "test_ovf32", dataformat = "binary")

JuMag.cuda_using_double(true)
sim = Sim(mesh, name = "test_ovf3")
set_Ms(sim, 8.0e5)
read_ovf(sim, "test_ovf")
println(sim.spin)

spin = Array(sim.spin)
@test isapprox(spin, expected,atol=1e-5)
spin = Array(sim.prespin)
@test isapprox(spin, expected,atol=1e-5)

read_ovf(sim, "test_ovf32")
println(sim.spin)
spin = Array(sim.spin)
@test isapprox(spin, expected, atol=1e-5)
spin = Array(sim.prespin)
@test isapprox(spin, expected, atol=1e-5)

JuMag.cuda_using_double(false)
sim = Sim(mesh, name = "test_ovf4")
set_Ms(sim, 8.0e5)
read_ovf(sim, "test_ovf")
println(sim.spin)
spin = Array(sim.spin)
@test isapprox(spin, expected,atol=1e-5)
spin = Array(sim.prespin)
@test isapprox(spin, expected,atol=1e-5)

relax(sim,maxsteps = 20)

read_ovf(sim, "test_ovf32")
println(sim.spin)
spin = Array(sim.spin)
@test isapprox(spin, expected, atol=1e-5)
spin = Array(sim.prespin)
@test isapprox(spin, expected, atol=1e-5)
