using SpinDynamics
using Test

#Test mesh
mesh =  FDMesh(dx=1.1, nx=10)
@test mesh.dx == 1.1
@test mesh.nx == 10
#println(mesh.ngbs)
@test mesh.ngbs[1,1] == -1
@test mesh.ngbs[1,10] == 9

sim = Sim(mesh, driver="LLG")

@test sim.nxyz == 10

init_m0(sim, (1.0, 1.0, 0))

@test isapprox(sim.spin[1],0.5*2^0.5)
@test isapprox(sim.spin[2],0.5*2^0.5)
@test sim.spin[3] == 0.0
println(sim.spin)
