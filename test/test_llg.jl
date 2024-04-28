using JuMag
using Test

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_llg(integrator="DormandPrince")
    #Test mesh
    mesh = FDMesh(; nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.tol = 1e-8

    add_zeeman(sim, (0, 0, 1e5))

    init_m0(sim, (1.0, 0, 0))

    for i in 1:300
        run_until(sim, 1e-12 * i)
    end

    #println(sim.spin[1]," ",sim.spin[2]," ",sim.spin[3])
    ts = Array([3e-10])
    mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)

    m = Array(sim.spin)

    @test abs(mx[1] - m[1]) < 8e-7
    @test abs(my[1] - m[2]) < 8e-7
    @test abs(mz[1] - m[3]) < 8e-7
end


function test_llg_rk(integrator="RungeKutta")
    #Test mesh
    mesh = FDMesh(; nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)
    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.dt = 1e-12

    add_zeeman(sim, (0, 0, 1e5))
    init_m0(sim, (1.0, 0, 0))

    for i in 1:300
        JuMag.run_step(sim, sim.driver)
    end

    #println(sim.spin[1]," ",sim.spin[2]," ",sim.spin[3])
    ts = Array([3e-10])
    mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)

    m = Array(sim.spin)

    @test abs(mx[1] - m[1]) < 8e-7
    @test abs(my[1] - m[2]) < 8e-7
    @test abs(mz[1] - m[3]) < 8e-7
end

@testset "Test LLG CPU" begin
    set_backend("cpu")
    for integrator in ["DormandPrince", "DormandPrinceCayley"]
        test_llg(integrator)
    end

    for integrator in ["RungeKutta", "RungeKuttaCayley"]
        test_llg_rk(integrator)
    end    
end

@testset "Test LLG CUDA" begin
    if Base.find_package("CUDA") !== nothing
        using CUDA
        for integrator in ["DormandPrince", "DormandPrinceCayley"]
            test_llg(integrator)
        end
    end
end
