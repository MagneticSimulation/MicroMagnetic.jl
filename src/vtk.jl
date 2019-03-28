using WriteVTK

function save_vtk(sim::AbstractSim, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  xyz = zeros(Float32, 3, nx+1, ny+1, nz+1)
  dx, dy, dz = 1.0,1.0,1.0
  if isa(sim, MicroSim)
    dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
  else
    dx, dy, dz=mesh.a, mesh.a, mesh.a
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
        b = reshape(sim.field, (3, nx, ny, nz))
        vtk_cell_data(vtk, b, i.name)
      end
    end
  end
  vtk_save(vtk)
end
