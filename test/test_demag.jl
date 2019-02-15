using SpinDynamics
using Test

mesh =  create_mesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=2, ny=1, nz=1)
Ms = 8.6e5
sim = create_sim(mesh)
sim.Ms[:] .= Ms

init_m0(sim, (1,0,0))
add_demag(sim)

SpinDynamics.effective_field(sim, sim.spin, 0.0)
println(sim.field)
@test isapprox(sim.field, [-170551.8913984, 0,0,-170551.8913984,0,0])



mesh =  create_mesh(nx=3, ny=2, nz=1)
Ms = 8.6e5
sim = create_sim(mesh)
sim.Ms[:] .= Ms
init_m0(sim, (0.1,0.2,1))
add_demag(sim)

SpinDynamics.effective_field(sim, sim.spin, 0.0)
println(sim.field)
expected = [-9615.99019074, -39898.78767025, -430282.70478141,-8664.33293854,
  -51323.59117349,-496012.77748287, -27749.99302205,-48965.78908591,
 -430282.70478141,-27749.99302205,  -48965.78908591, -430282.70478141,
   -8664.33293854,  -51323.59117349, -496012.77748287,   -9615.99019074,
  -39898.78767025, -430282.70478141];
@test isapprox(sim.field, expected)
