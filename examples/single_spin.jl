using SpinDynamics
using Base.Test

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha*alpha)
    beta = precession * H0 * ts

    mx = cos(beta) ./ cosh(alpha .* beta)
    my = sin(beta) ./ cosh(alpha .* beta)
    mz = tanh(alpha .* beta)
    return mx, my, mz
end

#Test mesh
mesh =  SpinDynamics.create_mesh(nx=1)

sim = SpinDynamics.create_sim(mesh, name="dyn")

sim.alpha = 0.5
#sim.mu_s = 1.0
sim.Hz = 1.0

ts = [0.2*i for i in 1:10]
mx, my = SpinDynamics.solver(sim, ts)
print(mx)

#using PlotlyJS
using PyPlot
ts = [0.2*i for i in 0:10]
ay1, ay2, ay3 = analytical(0.5, 1.0, 1.0, ts)
plot(ts, mx, color=:red, linewidth=2, label="sim")
plot(ts, ay1, color=:blue, linewidth=2, label="al", dashes=(2,2))
savefig("m_ts.png")
