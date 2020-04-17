abstract type Geometry end

@enum Axis begin
          ex = 1
          ey = 2
          ez = 3
      end

mutable struct Cylinder <: Geometry
    axis::Axis
    r_a::Float64
    r_b::Float64
    h::Float64
    xc :: Float64
    yc :: Float64
    zc :: Float64
    shape :: Array{Bool}
    Cylinder() = new()
end

function create_cylinder(mesh::Mesh, axis::Axis)
    geo = Cylinder()
    geo.axis = axis
    geo.xc =  mesh.nx*mesh.dx/2.0  #[0, nx*dx]
    geo.yc =  mesh.ny*mesh.dy/2.0  #[0, ny*dy]
    geo.zc =  mesh.nz*mesh.dz/2.0  #[0, nz*dz]
    geo.shape = zeros(Bool, mesh.nxyz)
    if axis == ez
        r_a = geo.xc
        r_b = geo.yc
        h = 2*geo.zc
        for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
            r1 = (i-0.5)*mesh.dx - geo.xc
            r2 = (j-0.5)*mesh.dy - geo.yc
            id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
            if (r1/r_a)^2+(r2/r_b)^2<=1
                geo.shape[id] = true
            end
        end
    elseif axis == ex
        r_a = geo.yc
        r_b = geo.zc
        h = 2*geo.xc
        for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
            r1 = (j-0.5)*mesh.dy - geo.yc
            r2 = (k-0.5)*mesh.dz - geo.zc
            id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
            if (r1/r_a)^2+(r2/r_b)^2<=1
                geo.shape[id] = true
            end
        end
    elseif axis == ey
        r_a = geo.zc
        r_b = geo.xc
        h = 2*geo.yc
        for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
            r1 = (k-0.5)*mesh.dz - geo.zc
            r2 = (i-0.5)*mesh.dx - geo.xc
            id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
            if (r1/r_a)^2+(r2/r_b)^2<=1
                geo.shape[id] = true
            end
        end
    end
    return geo
end
