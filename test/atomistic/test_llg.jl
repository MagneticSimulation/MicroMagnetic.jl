using NuMag
using Test

function analytical_llg(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha*alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_llg()
    #Test mesh
    mesh =  CubicMesh(nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh, name="spin")

    set_mu_s(sim, 2*mu_B)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5/mu_0
    sim.driver.integrator.tol = 1e-8

    add_zeeman(sim, (0, 0, 0.1)) #Tesla

    init_m0(sim, (1.0, 0, 0))

    ts = Float64[]
    mx = Float64[]
    my = Float64[]
    mz = Float64[]
    for i=1:300
      run_until(sim, 1e-12*i)
    end

    println(sim.driver.integrator.t)
    spin = Array(sim.spin)

    #println(sim.spin[1]," ",sim.spin[2]," ",sim.spin[3])
    ts = Array([3e-10])
    mx, my, mz = analytical_llg(0.05, 2.21e5/mu_0, 0.1, ts)
    println(mx[1]-spin[1])
    println(my[1]-spin[2])
    println(mz[1]-spin[3])

    @test abs(mx[1]-spin[1]) < 8e-7
    @test abs(my[1]-spin[2]) < 8e-7
    @test abs(mz[1]-spin[3]) < 8e-7
end


@testset "Test LLG CPU" begin
  set_backend("cpu")
  test_llg()
end

@testset "Test LLG CUDA" begin
  if Base.find_package("CUDA") !== nothing
      using CUDA
      test_llg()
  end
end
