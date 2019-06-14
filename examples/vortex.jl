using JuMag
using Printf
using NPZ

JuMag.using_gpu()
JuMag.cuda_using_double(true)
mesh =  FDMeshGPU(dx=2e-9, dy=2e-9, dz=5e-9, nx=100, ny=100, nz=4)
#mesh = FDMesh(dx=2e-9, dy=2e-9, dz=3e-9, nx=100, ny=100, nz=4)

function circular_Ms(i,j,k,dx,dy,dz)
    x = i-50.5
    y = j-50.5
    r = (x^2+y^2)^0.5
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8e5
    end
    return 0.0
end

function init_fun(i,j,k,dx,dy,dz)
  x = i-50.5
  y = j-50.5
  r = (x^2+y^2)^0.5
  if r<5
    return (0,0,1)
  end
  return (y/r, -x/r, 0)
end

function relax_system()
  sim = Sim(mesh, driver="SD", name="sim")
  set_Ms(sim, circular_Ms)

  add_exch(sim, 1.3e-11, name="exch")
  add_demag(sim)

  init_m0(sim, init_fun)
  relax(sim, maxsteps=2000, stopping_torque=1.0, save_vtk_every = 100, save_m_every=-1)
  #println(sim.spin)
  npzwrite("m0.npy", sim.spin)
  #save_vtk(sim, "skx", fields=["exch", "dmi"])
end

relax_system()
