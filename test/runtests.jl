using JuMag

include("test_mesh.jl")
include("test_zeeman.jl")
include("test_anis.jl")

@testset "CPU back-end" begin
  test_FDMesh()
  test_zeeman()
  test_anis()
  test_cubic_anis()
end

@testset "NVIDIA CUDA back-end" begin
  if Base.find_package("CUDA") !== nothing
    using CUDA
    test_FDMesh()
    test_zeeman()
    test_anis()
    test_cubic_anis()
  end
end


