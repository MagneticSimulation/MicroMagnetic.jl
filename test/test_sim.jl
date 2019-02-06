using SpinDynamics
using Test

#Test mesh
mesh =  SpinDynamics.create_mesh(dx=1.1, nx=10)
@test mesh.dx==1.1
@test mesh.nx==10
#println(mesh.ngbs)
@test mesh.ngbs[1,1] == -1
@test mesh.ngbs[1,10] == 9

sim = SpinDynamics.create_sim(mesh)
println(sim.nxyz)

sim.spin[1,1]=1.0
println(sim.spin)

SpinDynamics.compute_exch(sim.spin, sim.field, sim.energy, mesh.ngbs, 1.0)
println(sim.field)
