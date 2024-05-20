using MicroMag
using Test

MicroMag.cuda_using_double(true)

mesh = FDMeshGPU(dx=2e-9, dy=2e-9, dz=2e-9, nx=1, ny=1, nz=3, pbc="xy")

function m0_anti(i,j,k, dx, dy, dz)
  if k == 1
    return (0,0,1)
  end
  if k == 3
    return (0,0,-1)
  end
  return (1,1,1)
end

function spatial_Ms(i,j,k, dx, dy, dz)
  if k == 1 || k== 3
    return 8e5
  end
  return 0
end

sim = Sim(mesh)
set_Ms(sim, spatial_Ms)
init_m0(sim, m0_anti)
add_exch(sim, 1.3e-11)
#add_anis(sim, 5e4, axis=(0,0,1))
#add_exch_rkky(sim, -1e-4)
add_dmi_interlayer(sim, (0.0, 1e-4, 0))

relax(sim, stopping_dmdt=0.01)

expected = [-sqrt(2)/2, 0, sqrt(2)/2, 0,0,0, -sqrt(2)/2, 0, -sqrt(2)/2]
println(maximum(Array(sim.spin) .- expected))
@test maximum(Array(sim.spin) .- expected) < 2e-5


