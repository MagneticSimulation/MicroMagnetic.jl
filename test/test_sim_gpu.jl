using SpinDynamics
using Test

#Test mesh
mesh =  FDMeshGPU(dx=1.1e-9, nx=10)
@test mesh.dx == SpinDynamics.FloatGPU(1.1e-9)
@test mesh.nx == 10
#println(mesh.ngbs)
@test mesh.volume == SpinDynamics.FloatGPU(1.1e-27)
