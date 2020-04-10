using JuMag
using Test

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha*alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_integrator(;integrator="RungeKutta")
    mesh =  FDMesh(nx=1, ny=1, dx=1e-9)
    sim = Sim(mesh, name="spin", integrator=integrator)
    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    if integrator=="Default"
        sim.driver.ode.tol = 1e-6
    end
    add_zeeman(sim, (0, 0, 1e5))

    init_m0(sim, (1.0, 0, 0))

    for i=1:300
        JuMag.advance_step(sim, sim.driver.ode)
    end
    t = sim.driver.ode.t
    println("Running at time=", t)

    run_until(sim, t+5e-10)
    t = sim.driver.ode.t
    println("Running at time=", t)

    ts = Array([t])
    mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)

    @test abs(mx[1]-sim.spin[1]) < 8e-6
    @test abs(my[1]-sim.spin[2]) < 8e-6
    @test abs(mz[1]-sim.spin[3]) < 8e-6
end

test_integrator(integrator="Default")
println("Test 'Default' Done!")
#test_integrator(integrator="RungeKutta")
#println("Test 'RungeKutta' Done!")
test_integrator(integrator="DormandPrince")
println("Test 'DormandPrince' Done!")
