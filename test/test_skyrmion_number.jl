using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
  include("test_utils.jl")
end

function m0_fun(i,j,k,dx,dy,dz)
  r2 = (i-50)^2 + (j-50)^2
  if r2 < 10
    return (0.1, 0, -1)
  end
  return (0,0,1)
end

function test_skyrmion_number()
  mesh =  FDMesh(nx=100, ny=100, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")
  sim = Sim(mesh, driver="SD", name="sim")
  set_Ms(sim, 8e5)

  add_exch(sim, 1.3e-11, name="exch")
  add_zeeman(sim, (0,0,4e5))
  add_dmi(sim, 4e-3, name="dmi")

  init_m0(sim, m0_fun)
  relax(sim, maxsteps=2000, stopping_dmdt=0.1, save_m_every=-1)

  m = Array(sim.spin)
  v = zeros(eltype(m), sim.n_total)
  compute_skyrmion_number(v, m, mesh)
  Rxs, Rys = compute_guiding_center(m, mesh)
  println(sum(v)," ", Rxs, " ", Rxs)
  Q = sum(v)
  @test abs(Q+1) < 1e-6
  @test abs(Rxs[1]-100e-9)<1e-9
  @test abs(Rys[1]-100e-9)<1e-9

end

@using_gpu()
test_functions("skx_number", test_skyrmion_number)