using JuMag
using Printf
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
mesh =  FDMesh(nx=1, dx=1, unit_length=1e-9)

sim = Sim(mesh, name="spin")

set_Ms(sim, 8e5)
sim.driver.alpha = 0.1
sim.driver.gamma = 2.21e5

add_zeeman(sim, (0, 0, 1e5))

init_m0(sim, (1.0, 0, 0))

ts = Float64[]
mx = Float64[]
my = Float64[]
mz = Float64[]
for i=1:100
  run_until(sim, 1e-11*i)
  @info @sprintf "t=%g sim.t=%g mx=%g my=%g mz=%g" i*1e-11 sim.driver.ode.t sim.spin[1] sim.spin[2] sim.spin[3];
  push!(ts, i*1e-11)
  push!(mx, sim.spin[1])
  push!(my, sim.spin[2])
  push!(mz, sim.spin[3])
end
#JuMag.advance_step(sim, sim.ode)

#using PlotlyJS
using Plots
gr()

#ay1, ay2, ay3 = analytical(0.1, 2.21e5, 1e5, ts)
ay1, ay2, ay3 = analytical(0.1, 2.21e5, 1e5, ts)
plot(ts, mx, marker=:h, markersize=3, linewidth=1, label="sim")
plot!(ts, my, marker=:c, markersize=3, linewidth=1, label="sim")
plot!(ts, mz, marker=:p, markersize=1.2, linewidth=1, label="sim")
plot!(ts, ay1, color=:blue, label="mx")
plot!(ts, ay2, color=:blue, label="my")
plot!(ts, ay3, color=:blue, label="mz")
savefig("mxyz.png")
