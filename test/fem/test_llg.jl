using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha*alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_llg(;integrator="DormandPrince")

  mesh = FEMesh("meshes/one_hex.mesh", unit_length=1e-9)
  sim = Sim(mesh, name="spin", integrator=integrator)

  set_Ms(sim, 8e5)
  sim.driver.alpha = 0.05
  sim.driver.gamma = 2.21e5
  sim.driver.integrator.tol=1e-9

  add_zeeman(sim, (0, 0, 1e5))

  init_m0(sim, (1.0, 0, 0))

  t_stop = 1e-9
  run_until(sim, t_stop)

  spin =  Array(sim.spin)

  println(spin[1]," ",spin[2]," ",spin[3])
  ts = Array([t_stop])
  mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)
  println(mx[1] - spin[1])
  println(my[1] - spin[2])
  println(mz[1] - spin[3])

  @test abs(mx[1] - spin[1]) < 8e-7
  @test abs(my[1] - spin[2]) < 8e-7
  @test abs(mz[1] - spin[3]) < 8e-7
end


@using_gpu()
test_functions("LLG (FE)", test_llg)