using MicroMagnetic
using Test

include("test_mesh.jl")
include("test_shapes.jl")
include("test_zeeman.jl")
include("test_anis.jl")
include("test_exch.jl")
include("test_dmi.jl")
include("test_demag.jl")
include("test_fields.jl")
include("test_sim.jl")
include("test_llg.jl")
include("test_relax.jl")
include("test_pins.jl")
include("test_stt.jl")
include("test_ovf.jl")
include("test_interlayer.jl")
include("test_skyrmion_number.jl")

include("atomistic/test_mesh.jl")
include("atomistic/test_llg.jl")
include("atomistic/test_fields.jl")

include("test_mc.jl")
