using JuMag

include("test_mesh.jl")

@testset "CPU back-end" begin
  test_FDMesh()
end



@testset "NVIDIA CUDA back-end" begin
  if Base.find_package("CUDA") !== nothing
    using CUDA

    test_FDMesh()
  end
end


