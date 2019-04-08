using JuMag
using Test

mesh =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=2, ny=1, nz=1)
Ms = 8.6e5
sim = Sim(mesh)

set_Ms(sim, Ms)

init_m0(sim, (1,0,0))
add_demag(sim)

JuMag.effective_field(sim, sim.spin, 0.0)
expected = [-170551.8913984, 0,0,-170551.8913984,0,0]
@test isapprox(Array(sim.field), expected, rtol=1e-6)
println(sim.field)

mesh =  FDMeshGPU(nx=3, ny=2, nz=1)
Ms = 8.6e5
sim = Sim(mesh)
set_Ms(sim, Ms)
init_m0(sim, (0.1,0.2,1))
add_demag(sim)

JuMag.effective_field(sim, sim.spin, 0.0)

expected = [-9615.99019074, -39898.78767025, -430282.70478141,-8664.33293854,
  -51323.59117349,-496012.77748287, -27749.99302205,-48965.78908591,
 -430282.70478141,-27749.99302205,  -48965.78908591, -430282.70478141,
   -8664.33293854,  -51323.59117349, -496012.77748287,   -9615.99019074,
  -39898.78767025, -430282.70478141];
@test isapprox(Array(sim.field), expected, rtol=1e-6)
#println(Array(sim.field)-expected)

mesh =  FDMeshGPU(dx=2, dy=1, dz=1, nx=12, ny=1, nz=1)
Ms = 8.6e5
sim = Sim(mesh)
set_Ms(sim, Ms)

init_m0(sim, (0.5,0.6,0.7))
add_demag(sim)

JuMag.effective_field(sim, sim.spin, 0.0)
expected = [ -40715.4448514, -221564.08111452, -258491.42796693,   -3885.3038712,
 -243662.16570271, -284272.52665316 ,  -1420.95983327, -245140.77212539,
 -285997.56747963,   -785.62732964, -245521.97162757, -286442.30023216,
    -550.58082889, -245662.99952803, -286606.8327827,     -464.3653603,
 -245714.72880919, -286667.18361072,    -464.3653603,  -245714.72880919,
 -286667.18361072,    -550.58082889, -245662.99952803, -286606.8327827,
    -785.62732964, -245521.97162757, -286442.30023216,   -1420.95983327,
 -245140.77212539, -285997.56747963,   -3885.3038712,  -243662.16570271,
 -284272.52665316,  -40715.4448514,  -221564.08111452, -258491.42796693];
@test isapprox(Array(sim.field), expected, rtol=1e-6)
#println(Array(sim.field)-expected)
