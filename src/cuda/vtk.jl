using WriteVTK
using NPZ

function save_m(sim::AbstractSimGPU, fname::String; vtk::Bool = false, npy::Bool = false, vtk_folder::String="vtks",npy_folder::String="npys")
  if vtk
      !isdir(vtk_folder) && mkdir(vtk_folder)
      save_vtk(sim, joinpath(vtk_folder, fname))
  end
  if npy
    !isdir(npy_folder) && mkdir(npy_folder)
    save_npy(sim,joinpath(npy_folder,fname))
  end
end

function save_npy(sim::AbstractSimGPU,name::String)
  T = _cuda_using_double.x ? Float64 : Float32
  name = @sprintf("%s.npy", name)
  mesh = sim.mesh
  n_total = mesh.nx*mesh.ny*mesh.nz
  spin = zeros(T, 3*n_total)
  copyto!(spin, sim.spin)
  npzwrite(name, spin)
end

function save_vtk_triangular_mesh(sim::AbstractSimGPU, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  xyz = zeros(Float32, 3, nx, ny, nz)
  dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
  for k = 1:nz, j = 1:ny, i = 1:nx
    xyz[1, i, j, k] = (i-0.5)*dx + (j-1)*dx/2
    xyz[2, i, j, k] = (j-0.5)*dy
    xyz[3, i, j, k] = (k-0.5)*dz
  end

  vtk = vtk_grid(fname, xyz)
  T = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(T, 3*sim.n_total)
  copyto!(spin, sim.spin)
  b = reshape(spin, (3, nx, ny, nz))
  vtk_point_data(vtk, b , "m")

  if length(fields) > 0
    compute_fields_to_gpu(sim,sim.spin,0.0)
    fields = Set(fields)
    for i in sim.interactions
      if i.name in fields
        b = reshape(i.field, (3, nx, ny, nz))
        vtk_point_data(vtk, b, i.name)
      end
    end
  end
  vtk_save(vtk)
end

function save_vtk_fcc_mesh(sim::AbstractSimGPU, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh

  vtk = vtk_grid(fname, mesh.coordinates, MeshCell[])

  T = _cuda_using_double.x ? Float64 : Float32
  vtk_point_data(vtk, Array(sim.spin), "m")

  if length(fields) > 0
    compute_fields_to_gpu(sim,sim.spin,0.0)
    fields = Set(fields)
    for i in sim.interactions
      if i.name in fields
        vtk_point_data(vtk, Array(i.field), i.name)
      end
    end
  end
  vtk_save(vtk)
end


function save_vtk(sim::AbstractSimGPU, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh

  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  xyz = zeros(Float32, 3, nx+1, ny+1, nz+1)
  dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
  
  if isa(mesh, TriangularMeshGPU)
      save_vtk_triangular_mesh(sim, fname, fields=fields)
      return nothing
  elseif isa(mesh, FccMeshGPU)
      save_vtk_fcc_mesh(sim, fname, fields=fields)
      return nothing
  else
    for k = 1:nz+1, j = 1:ny+1, i = 1:nx+1
      xyz[1, i, j, k] = (i-0.5-nx/2)*dx
      xyz[2, i, j, k] = (j-0.5-ny/2)*dy
      xyz[3, i, j, k] = (k-0.5-nz/2)*dz
    end
  end
  vtk = vtk_grid(fname, xyz)
  T = _cuda_using_double.x ? Float64 : Float32
  xyz = zeros(T, 3*sim.n_total)
  copyto!(xyz, sim.spin)
  b = reshape(xyz, (3, nx, ny, nz))
  vtk_cell_data(vtk, b , "m")

  if length(fields) > 0
    compute_fields_to_gpu(sim,sim.spin,0.0)
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



function save_vtk_points(sim::AbstractSimGPU, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh

  if isa(mesh, FccMeshGPU)
    save_vtk_fcc_mesh(sim, fname, fields=fields)
    return nothing
  end

  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  xyz = zeros(Float32, 3, nx, ny, nz)
  dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
  for k = 1:nz, j = 1:ny, i = 1:nx
    xyz[1, i, j, k] = (i-0.5-nx/2)*dx
    xyz[2, i, j, k] = (j-0.5-ny/2)*dy
    xyz[3, i, j, k] = (k-0.5-nz/2)*dz
  end
  vtk = vtk_grid(fname, xyz)
  T = _cuda_using_double.x ? Float64 : Float32

  m = zeros(T, 3*sim.n_total)
  copyto!(m, sim.spin)

  b = reshape(m, (3, nx, ny, nz))
  vtk_point_data(vtk, b , "m")

  if length(fields) > 0
    compute_fields_to_gpu(sim,sim.spin,0.0)
    fields = Set(fields)
    for i in sim.interactions
      if i.name in fields
        b = reshape(i.field, (3, nx, ny, nz))
        vtk_point_data(vtk, b, i.name)
      end
    end
  end
  vtk_save(vtk)
end

# only works for CylindricalTubeMeshGPU
function save_vtu(sim::AbstractSimGPU, fname::String; fields::Array{String, 1} = String[])
  mesh = sim.mesh
  points = mesh.coordinates
  ngbs = Array(mesh.ngbs)
  cells = MeshCell[]
  for k = 1:mesh.nz-1, i=1:mesh.nr
      id = index(i, 1, k, mesh.nr, 1, mesh.nz)
      id2 = indexpbc(i+1, 1, k+1, mesh.nr, 1, mesh.nz, true, false, false)
      cell = MeshCell(VTKCellTypes.VTK_QUAD, [id, ngbs[2, id], id2, ngbs[4, id]])
      push!(cells, cell)
  end

  vtk_grid(fname, points, cells) do vtk
      m =  Array(sim.spin)
      mc = convert_m_to_cylindrical(m, mesh.nr, mesh.nz)
      vtk_point_data(vtk, m, "m"; component_names=["mx", "my", "mz"])
      vtk_point_data(vtk, mc, "mc"; component_names=["mr", "mt", "mz"])

      if length(fields) > 0
        compute_fields_to_gpu(sim,sim.spin,0.0)
        fields = Set(fields)
        for i in sim.interactions
          if i.name in fields
            vtk_point_data(vtk, i.field, i.name)
          end
        end
      end

  end
end


