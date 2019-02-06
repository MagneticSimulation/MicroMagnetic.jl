using SpinDynamics
using Test

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha*alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

#Test mesh
mesh =  SpinDynamics.create_mesh(nx=1)

sim = SpinDynamics.create_sim(mesh, name="dyn")

sim.alpha = 0.1
sim.gamma = 2.21e5
#sim.Ms = 1.0
sim.Hz = 1e5

SpinDynamics.init_m0(sim, 1.0, 0.0, 0.0)

ts = Float64[]
mx = Float64[]
my = Float64[]
mz = Float64[]
for i=1:100
  SpinDynamics.run_until(sim, 1e-11*i)
  println(i, sim.ode.t, sim.spin)
  push!(ts, i*1e-11)
  push!(mx, sim.spin[1])
  push!(my, sim.spin[2])
  push!(mz, sim.spin[3])
end
#SpinDynamics.advance_step(sim, sim.ode)

#using PlotlyJS
using Plots
gr()

ay1, ay2, ay3 = analytical(0.1, 2.21e5, 1e5, ts)
plot!(ts, mx, line=(:dot, 2), marker=([:hex :d]), label="sim")
plot!(ts, my, line=(:dot, 2), markersize=2, label="sim")
plot!(ts, mz, line=(:dot, 2), label="sim")
plot!(ts, ay1, color=:blue, label="al")
plot!(ts, ay1, color=:blue, label="al")
plot!(ts, ay1, color=:blue, label="al")
savefig("m_ts.png")
