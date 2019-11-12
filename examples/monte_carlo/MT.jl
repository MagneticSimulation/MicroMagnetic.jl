using JuMag
using Random

function rand_m(i, j, k)
  return 2*rand(3).-1
end

function rand_m_uniform(i, j, k)
  return (0,0,1)
end

JuMag.cuda_using_double(true)

function relax_system(T)
  mesh =  CubicMeshGPU(nx=30, ny=30, nz=30, pbc="xyz")
  sim = MonteCarloNew(mesh, name="mc")
  init_m0(sim, rand_m_uniform)

  add_exch(sim, J=300*k_B)
  add_dmi(sim, D=0, D1=0)
  add_zeeman(sim, Hx=0, Hy=0, Hz=0)
  add_anis(sim, Ku=0, Kc=0)

  sim.T = 100000
  run_sim(sim, maxsteps=50000, save_vtk_every=-1, save_m_every=-1)
  sim.T = T
  run_sim(sim, maxsteps=200000, save_vtk_every=-1, save_m_every=-1)

  ms = zeros(1000)
  sim.T = T
  for i = 1:1000
      run_sim(sim, maxsteps=100, save_vtk_every=-1, save_m_every=-1)
      t = JuMag.average_m(sim)
      ms[i] = sqrt(t[1]^2+t[2]^2+t[3]^2)
  end
  return sum(ms)/length(ms)
end

f = open("M_H.txt", "w")
write(f, "#T(K)     m \n")
for T = 10:20:500
    println("Running for $T ...")
    m = relax_system(T)
    write(f, "$T    $m \n")
end
close(f)
