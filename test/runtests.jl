using JuMag

include("test_mesh.jl")
include("test_zeeman.jl")

@testset "CPU back-end" begin
  test_FDMesh()
  test_zeeman()
end



@testset "NVIDIA CUDA back-end" begin
  if Base.find_package("CUDA") !== nothing
    using CUDA
    test_FDMesh()
    test_zeeman()
  end
end


