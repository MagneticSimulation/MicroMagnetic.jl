using JuMag

include("test_mesh.jl")
include("test_exch.jl")
include("test_demag.jl")
include("test_demag_pbc.jl")
include("test_zeeman.jl")
include("test_interface.jl")
include("test_fields.jl")
include("test_stt.jl")
include("test_llg.jl")
include("test_relax.jl")
include("test_util.jl")
include("test_ovf.jl")

if JuMag._cuda_available.x
  JuMag.cuda_using_double()
  include("cuda/test_sim.jl")
  include("cuda/test_exch.jl")
  include("cuda/test_demag.jl")
  include("cuda/test_llg.jl")
  include("cuda/test_ovf.jl")

  mesh =  FDMeshGPU(nx=500, ny=1, nz=11, dx=2e-9, dy=2e-9, dz=1e-9)
  relax_system(mesh)
  for (beta, u) in [(0, 10), (0.1, 3.2), (0.2, 4.7)]
    run_dynamics_stt(mesh, beta=beta, u=u)
  end
  relax_system_SD(mesh)
  relax_system_LLG(mesh)

  JuMag.cuda_using_double()
  test_zeeman(gpu=true)
  mesh =  FDMeshGPU(nx=20, ny=5, nz=3, dx=2.5e-9, dy=2.5e-9, dz=3e-9)
  test_fields(mesh, gpu=true)
end
