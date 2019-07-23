using JuMag
using Test

mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=1)
Ms = 8.6e5
sim = Sim(mesh)
sim.Ms[:] .= Ms

init_m0(sim, (0,0,1))
add_demag(sim, Nx=20, Ny=20)

JuMag.effective_field(sim, sim.spin, 0.0)
println(sim.field)
#@test isapprox(sim.field, [-170551.8913984, 0,0,-170551.8913984,0,0])
