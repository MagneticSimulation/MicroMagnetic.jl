using JuMag
using Test

#Test mesh
mesh =  FDMeshGPU(dx=1.1e-9, nx=10, ny=2)
@test isapprox(mesh.dx*1.0, 1.1e-9, rtol=1e-7) #FIXME: shall we always use Float64 for dx, dy and dz?
@test mesh.nx == 10
#println(mesh.ngbs)
@test isapprox(mesh.volume*1.0, 1.1e-27, rtol=1e-7)

sim = Sim(mesh, driver="LLG")
set_Ms(sim, 1.0)

@test sim.nxyz == 20

init_m0(sim, (1.0, 1.0, 0))
println(sim.spin)
spin = Array(sim.spin)

@test isapprox(spin[1],0.5*2^0.5)
@test isapprox(spin[2],0.5*2^0.5)


set_Ms(sim, 8.6e5)
sim.driver.alpha = 0.5
sim.driver.gamma = 2.21e5
sim.driver.precession = false

add_exch(sim, 1.3e-11)
JuMag.effective_field(sim, sim.spin, 0.0)
println(sim.spin)
