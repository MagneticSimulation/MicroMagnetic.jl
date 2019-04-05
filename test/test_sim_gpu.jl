using SpinDynamics
using Test

#Test mesh
mesh =  FDMeshGPU(dx=1.1e-9, nx=10, ny=2)
@test mesh.dx == SpinDynamics.FloatGPU(1.1e-9)
@test mesh.nx == 10
#println(mesh.ngbs)
@test mesh.volume == SpinDynamics.FloatGPU(1.1e-27)

sim = Sim(mesh, driver="LLG")

@test sim.nxyz == 20

init_m0(sim, (1.0, 1.0, 0))

spin = Array(sim.spin)

@test isapprox(spin[1],0.5*2^0.5)
@test isapprox(spin[2],0.5*2^0.5)


set_Ms(sim, 8.6e5)
sim.driver.alpha = 0.5
sim.driver.gamma = 2.21e5
sim.driver.precession = false

add_exch(sim, 1.3e-11)
SpinDynamics.effective_field(sim, sim.spin, 0.0)
