using JuMag

include("test_mesh.jl")
include("test_exch.jl")
include("test_demag.jl")
include("test_demag_pbc.jl")
include("test_zeeman.jl")
include("test_anis.jl")
include("test_fields.jl")
include("test_stt.jl")
include("test_llg.jl")
include("test_integrator.jl")
include("test_relax.jl")
include("test_util.jl")
include("test_ovf.jl")
include("test_neb.jl")

if JuMag._cuda_available.x
  JuMag.cuda_using_double()
  include("cuda/test_sim.jl")
  include("cuda/test_exch.jl")
  include("cuda/test_demag.jl")
  include("cuda/test_llg.jl")
  include("cuda/test_ovf.jl")
  include("cuda/test_neb.jl")

  include("atomistic/test_llg.jl")

  JuMag.cuda_using_double()
  test_zeeman(gpu=true)
  mesh =  FDMeshGPU(nx=20, ny=5, nz=3, dx=2.5e-9, dy=2.5e-9, dz=3e-9)
  test_fields(mesh, gpu=true)
end
