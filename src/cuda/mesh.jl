struct FDMeshGPU{T <: AbstractFloat} <: MeshGPU
  dx::T
  dy::T
  dz::T
  nx::Int64
  ny::Int64
  nz::Int64
  n_total::Int64
  xperiodic::Bool
  yperiodic::Bool
  zperiodic::Bool
  volume::T
end

"""
    FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1, pbc="open")

The GPU version of the FDMesh, which is needed for GPU compuation.
"""
function FDMeshGPU(;dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1, pbc="open")
  n_total = nx*ny*nz
  volume = dx*dy*dz
  Float = _cuda_using_double.x ? Float64 : Float32
  xperiodic = 'x' in pbc ? true : false
  yperiodic = 'y' in pbc ? true : false
  zperiodic = 'z' in pbc ? true : false
  return FDMeshGPU(Float(dx), Float(dy), Float(dz), nx, ny, nz, n_total, xperiodic, yperiodic, zperiodic, Float(volume))
end
