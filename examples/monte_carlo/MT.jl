using JuMag
using Random

function rand_m(i, j, k, nx, ny, nz)
  return 2*rand(3).-1
end

function rand_m_uniform(i, j, k, nx, ny, nz)
  return (0,0,1)
end

JuMag.cuda_using_double(true)

function relax_system()
  mesh =  CubicMeshGPU(nx=30, ny=30, nz=30, pbc="xyz")
  sim = MonteCarlo(mesh, name="mc")
  init_m0(sim, rand_m_uniform)

  add_exch(sim, J=300*k_B)
  add_dmi(sim, D=0, D1=0)
  add_zeeman(sim, Hx=0, Hy=0, Hz=0)
  add_anis(sim, Ku=0, Kc=0)

  sim.T = 100000
  run_sim(sim, maxsteps=50000, save_vtk_every=-1, save_m_every=-1)

  ms = zeros(1000)

  f = open("M_H.txt", "w")
  write(f, "#T(K)     m \n")

  for T = 500:-20:10
      sim.T = T
      println("Running for $T ...")
      run_sim(sim, maxsteps=100000, save_vtk_every=-1, save_m_every=1000)

      for i = 1:1000
          run_sim(sim, maxsteps=20, save_vtk_every=-1, save_m_every=-1)
          t = JuMag.average_m(sim)
          ms[i] = sqrt(t[1]^2+t[2]^2+t[3]^2)
      end

      m = sum(ms)/length(ms)
      write(f, "$T    $m \n")
  end
  close(f)

  return
end

relax_system()
