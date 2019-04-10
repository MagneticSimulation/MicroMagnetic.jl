using JuMag

include("test_sim.jl")
include("test_exch.jl")
include("test_demag.jl")
include("test_interface.jl")

if JuMag._cuda_available.x
  include("cuda/test_sim.jl")
  include("cuda/test_exch.jl")
  include("cuda/test_demag.jl")
  include("cuda/test_llg.jl")
end
