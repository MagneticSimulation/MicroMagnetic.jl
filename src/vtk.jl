using WriteVTK
using NPZ

"""
    save_m(sim::AbstractSim, fname::String; vtk::Bool = false, npy::Bool = false, vtk_folder = "vtks")

    If vtk = true, save spins to dir ./vtks/ in vtk format;
    If npy = true, save spins to current directory in npy format;

    Example:

    ```julia
        save_m(sim, sim.name, vtk = true, npy= true)
    ```
"""
function save_m(sim::AbstractSim, fname::String; vtk::Bool=false, npy::Bool=false, vtk_folder::String="vtks", npy_folder::String="npys")
  if vtk
    !isdir(vtk_folder) && mkdir(vtk_folder)
    save_vtk(sim, joinpath(vtk_folder, fname))
  end

  if npy
    !isdir(npy_folder) && mkdir(npy_folder)
    name = @sprintf("%s.npy", fname)
    npzwrite(name, joinpath(npy_folder,sim.spin))
  end
end

"""
    save_vtk(sim::AbstractSim, fname::String; fields::Array{String, 1} = String[])

Save magnetization or other fields to vtk.

```julia
    save_vtk(sim, "m")
```
"""
function save_vtk(sim::AbstractSim, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  xyz = zeros(Float32, 3, nx+1, ny+1, nz+1)
  dx, dy, dz = 1.0,1.0,1.0
  if isa(mesh, CubicMesh)
    dx, dy, dz=mesh.a, mesh.a, mesh.a
  else
    dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
  end
  for k = 1:nz+1, j = 1:ny+1, i = 1:nx+1
    xyz[1, i, j, k] = (i-0.5)*dx
    xyz[2, i, j, k] = (j-0.5)*dy
    xyz[3, i, j, k] = (k-0.5)*dz
  end
  vtk = vtk_grid(fname, xyz)
  b = reshape(sim.spin, (3, nx, ny, nz))
  vtk_cell_data(vtk, b , "m")

  if length(fields) > 0
    fields = Set(fields)
    for i in sim.interactions
      if i.name in fields
        b = reshape(i.field, (3, nx, ny, nz))
        vtk_cell_data(vtk, b, i.name)
      end
    end
  end
  vtk_save(vtk)
end


function ovf2vtk(mesh, ovf_name, vtk_name)
    sim = Sim(mesh, driver="SD")
    read_ovf(sim, ovf_name)
    save_vtk(sim, vtk_name)
    return nothing
end
