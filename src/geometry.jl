abstract type Geometry end

@enum Axis begin
          ex = 1
          ey = 2
          ez = 3
      end

mutable struct Cylinder <: Geometry
    axis::Axis
    r1::Float64  #radius of one surface (+z, +y, +x)
    r2::Float64  #radius of the other surface (-z, -y, -x)
    h::Float64
    xc :: Float64
    yc :: Float64
    zc :: Float64
    shape :: Array{Bool}
    Cylinder() = new()
end

function create_cylinder(mesh::Mesh, axis::Axis; r1=0, r2=0)
    geo = Cylinder()
    geo.axis = axis
    geo.xc =  mesh.nx*mesh.dx/2.0  #[0, nx*dx]
    geo.yc =  mesh.ny*mesh.dy/2.0  #[0, ny*dy]
    geo.zc =  mesh.nz*mesh.dz/2.0  #[0, nz*dz]
    geo.shape = zeros(Bool, mesh.nxyz)
    if axis == ez
        R = min(geo.xc, geo.yc)
        geo.r1 = r1 > 0 ? r1 : R
        geo.r2 = r2 > 0 ? r2 : R
        h = 2*geo.zc
        slope = (geo.r1 - geo.r2)/h
        b = (geo.r1 + geo.r2)/2
        for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
            x = (i-0.5)*mesh.dx - geo.xc
            y = (j-0.5)*mesh.dy - geo.yc
            z = (k-0.5)*mesh.dz - geo.zc
            r = sqrt(x^2+y^2)
            id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
            if slope*z + b > r
                geo.shape[id] = true
            end
        end
    elseif axis == ex
        R = min(geo.yc, geo.zc)
        geo.r1 = r1 > 0 ? r1 : R
        geo.r2 = r2 > 0 ? r2 : R
        h = 2*geo.xc
        slope = (geo.r1 - geo.r2)/h
        b = (geo.r1 + geo.r2)/2
        for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
            x = (i-0.5)*mesh.dx - geo.xc
            y = (j-0.5)*mesh.dy - geo.yc
            z = (k-0.5)*mesh.dz - geo.zc
            r = sqrt(z^2+y^2)
            id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
            if slope*x + b > r
                geo.shape[id] = true
            end
        end
    elseif axis == ey
        R = min(geo.xc, geo.yc)
        geo.r1 = r1 > 0 ? r1 : R
        geo.r2 = r2 > 0 ? r2 : R
        h = 2*geo.yc
        slope = (geo.r1 - geo.r2)/h
        b = (geo.r1 + geo.r2)/2
        for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
            x = (i-0.5)*mesh.dx - geo.xc
            y = (j-0.5)*mesh.dy - geo.yc
            z = (k-0.5)*mesh.dz - geo.zc
            r = sqrt(x^2+z^2)
            id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
            if slope*y + b > r
                geo.shape[id] = true
            end
        end
    end
    return geo
end
