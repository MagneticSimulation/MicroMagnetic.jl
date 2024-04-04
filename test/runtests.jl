using JuMag

include("test_mesh.jl")
include("test_zeeman.jl")
include("test_anis.jl")
include("test_exch.jl")

@testset "CPU back-end" begin
  test_FDMesh()
  test_zeeman()
  test_anis()
  test_cubic_anis()
  test_exch_scalar()
end

@testset "NVIDIA CUDA back-end" begin
  if Base.find_package("CUDA") !== nothing
    using CUDA
    test_FDMesh()
    test_zeeman()
    test_anis()
    test_cubic_anis()
    test_exch_scalar()
  end
end


