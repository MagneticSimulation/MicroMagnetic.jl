using JuMag

include("test_sim.jl")
include("test_exch.jl")
include("test_demag.jl")
include("test_interface.jl")

if JuMag._cuda_available.x
  include("test_sim.jl")
  include("test_exch.jl")
  include("test_demag.jl")
  include("test_llg.jl")
end
