using JuMag

include("test_mesh.jl")
include("test_zeeman.jl")
include("test_anis.jl")
include("test_exch.jl")
include("test_dmi.jl")

@testset "CPU back-end" begin
  test_FDMesh()
  test_zeeman()
  test_anis()
  test_cubic_anis()
  test_exch_scalar()
  test_exch_vectors()
  test_bulk_dmi()
  test_interfacial_dmi()
end

@testset "NVIDIA CUDA back-end" begin
  if Base.find_package("CUDA") !== nothing
    using CUDA
    test_FDMesh()
    test_zeeman()
    test_anis()
    test_cubic_anis()
    test_exch_scalar()
    test_exch_vectors()
    test_bulk_dmi()
    test_interfacial_dmi()
  end
end


