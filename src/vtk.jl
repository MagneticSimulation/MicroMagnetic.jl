using WriteVTK
using ReadVTK
using Printf

export save_vtk, ovf2vtk, read_vtk

#TODO: tidy up save_vtk

"""
    save_m(sim::AbstractSim, fname::String; vtk::Bool=false, vtk_folder::String="vtks")

    If vtk = true, save spins to dir ./vtks/ in vtk format;

    Example:

    ```julia
        save_m(sim, sim.name, vtk = true)
    ```
"""
function save_m(sim::AbstractSim, fname::String; vtk::Bool=false, vtk_folder::String="vtks")
    if vtk
        !isdir(vtk_folder) && mkdir(vtk_folder)
        save_vtk(sim, joinpath(vtk_folder, fname))
    end
end

"""
    save_vtk(sim::AbstractSim, fname::String; fields::Array{String, 1} = String[])

Save magnetization or other fields to vtk.

```julia
    save_vtk(sim, "m")
```
"""
function save_vtk(sim::AbstractSim, fname::String; fields::Array{String,1}=String[])
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    xyz = zeros(Float32, 3, nx + 1, ny + 1, nz + 1)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    for k in 1:(nz + 1), j in 1:(ny + 1), i in 1:(nx + 1)
        xyz[1, i, j, k] = (i - 0.5 - nx / 2) * dx
        xyz[2, i, j, k] = (j - 0.5 - ny / 2) * dy
        xyz[3, i, j, k] = (k - 0.5 - nz / 2) * dz
    end
    vtk = vtk_grid(fname, xyz)
    spin = zeros(eltype(sim.spin), length(sim.spin))
    copyto!(spin, sim.spin)

    b = reshape(spin, (3, nx, ny, nz))
    vtk_cell_data(vtk, b, "m")

    if length(fields) > 0
        fields = Set(fields)
        for i in sim.interactions
            if i.name in fields
                MicroMagnetic.effective_field(i, sim, sim.spin, 0.0)
                f = isa(i.field, Array) ? i.field : Array(i.field)
                b = reshape(f, (3, nx, ny, nz))
                vtk_cell_data(vtk, b, i.name)
            end
        end
    end
    return vtk_save(vtk)
end

function save_vtk_points(sim::AbstractSim, fname::String; fields::Array{String,1}=String[])
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    xyz = zeros(Float32, 3, nx, ny, nz)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    for k in 1:nz, j in 1:ny, i in 1:nx
        xyz[1, i, j, k] = (i - 0.5 - nx / 2) * dx
        xyz[2, i, j, k] = (j - 0.5 - ny / 2) * dy
        xyz[3, i, j, k] = (k - 0.5 - nz / 2) * dz
    end
    vtk = vtk_grid(fname, xyz)

    m = isa(sim.spin, Array) ? sim.spin : Array(sim.spin)

    b = reshape(m, (3, nx, ny, nz))
    vtk_point_data(vtk, b, "m")

    if length(fields) > 0
        fields = Set(fields)
        for i in sim.interactions
            if i.name in fields
                MicroMagnetic.effective_field(i, sim, sim.spin, 0.0)
                f = isa(i.field, Array) ? i.field : Array(i.field)
                b = reshape(f, (3, nx, ny, nz))
                vtk_point_data(vtk, b, i.name)
            end
        end
    end
    return vtk_save(vtk)
end

function ovf2_skyrmion_number_vtk(ovfname::String, fname::String)
    ovf = read_ovf(ovfname)
    m = ovf.data
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    dx, dy, dz = ovf.xstepsize, ovf.ystepsize, ovf.zstepsize

    mesh = FDMesh(; nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)
    v = zeros(nx * ny * nz)
    compute_skyrmion_number(v, m, mesh)
    xyz = zeros(Float32, 3, nx, ny, nz)

    for k in 1:nz, j in 1:ny, i in 1:nx
        xyz[1, i, j, k] = (i - 0.5 - nx / 2) * dx
        xyz[2, i, j, k] = (j - 0.5 - ny / 2) * dy
        xyz[3, i, j, k] = (k - 0.5 - nz / 2) * dz
    end
    vtk = vtk_grid(fname, xyz)
    b = reshape(v, (nx, ny, nz))
    vtk_point_data(vtk, b, "skxNumber")
    return vtk_save(vtk)
end
"""
    ovf2vtk(ovf_name, vtk_name=nothing; point_data=false, box=noting)

Convert ovf file to vtk format. The data will be saved to points if point_data == true otherwise
the data will be saved to cells.

If box is not nothing, it should be a tuple. For instance, box = (nx1, nx2, ny1, ny2, nz1, nz2).
In this case, the generated vtk only contains the spins inside the box (including the boundary).

```julia
    ovf2vtk("my.ovf", "test.vts")
    ovf2vtk("my.ovf", point_data=true)
```
"""
function ovf2vtk(ovf_name, vtk_name=nothing; point_data=false, box=nothing)
    if vtk_name == nothing
        vtk_name = endswith(ovf_name, ".ovf") ? ovf_name[1:(end - 4)] : ovf_name
    end
    ovf = read_ovf(ovf_name)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    if box != nothing
        nx = box[2] - box[1] + 1
        ny = box[4] - box[3] + 1
        nz = box[6] - box[5] + 1
    end

    if point_data
        xyz = zeros(Float32, 3, nx, ny, nz)
    else
        xyz = zeros(Float32, 3, nx + 1, ny + 1, nz + 1)
    end

    dim = size(xyz)
    for k in 1:dim[4], j in 1:dim[3], i in 1:dim[2]
        xyz[1, i, j, k] = (i - 0.5 - nx / 2) * ovf.xstepsize
        xyz[2, i, j, k] = (j - 0.5 - ny / 2) * ovf.ystepsize
        xyz[3, i, j, k] = (k - 0.5 - nz / 2) * ovf.zstepsize
    end

    vtk = vtk_grid(vtk_name, xyz)
    data = reshape(ovf.data, (3, ovf.xnodes, ovf.ynodes, ovf.znodes))
    b = data
    if box !== nothing
        b = data[:, box[1]:box[2], box[3]:box[4], box[5]:box[6]]
    end
    if point_data
        vtk_point_data(vtk, b, "m")
    else
        vtk_cell_data(vtk, b, "m")
    end
    vtk_save(vtk)
    return nothing
end

function read_vtk(filename; field="m", point_data=false)
    vtk = VTKFile(filename)
    data = point_data ? get_point_data(vtk) : get_cell_data(vtk)
    #point = get_points(vtk)
    f = copy(get_data(data[field]))
    return reshape(f, length(f))
end
