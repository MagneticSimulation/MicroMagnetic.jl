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

mutable struct Box <: Geometry
    x1::Float64
    y1::Float64
    z1::Float64
    x2 :: Float64
    y2 :: Float64
    z2 :: Float64
    shape :: Array{Bool}
    Box() = new()
end

function update_scalar_geometry(v::Array{T,1}, geo::Geometry, u::Number) where T<:AbstractFloat
    lv,lg = length(v),length(geo.shape)
    if lv != lg
        @error("Geometry doesn't fit array!")
    end
    for i = 1:lv
        if geo.shape[i]
            v[i] = u
        end
    end
end

function update_scalar_geometry(v::CuArray{T,1}, geo::Geometry, u::Number) where T<:AbstractFloat
    lv,lg = length(v),length(geo.shape)
    if lv != lg
        @error("Geometry doesn't fit array!")
    end
    for i = 1:lv
        if geo.shape[i]
            v[i] = u
        end
    end
end

"""
    create_cylinder(mesh::Mesh, axis::Axis; r1=0, r2=0)

Create a cylinder with axis(one of ex,ey or ez) and radius from r1 to r2. For example:
```julia
    mesh = FDMesh(nx=10,ny=10,nz=10,dx=1e-9,dy=1e-9,dz=1e-9)
    myCylinder = create_cylinder(mesh, ez, r1=10e-9, r2=10e-9)  
```
Usage maybe refers to function "create_box".
"""
function create_cylinder(mesh::Mesh, axis::Axis; r1=0, r2=0, xc="unvalued", yc="unvalued", zc="unvalued")
    geo = Cylinder()
    geo.axis = axis
    geo.xc =  xc == "unvalued" ? mesh.nx*mesh.dx/2.0 : xc  #[0, nx*dx]
    geo.yc =  yc == "unvalued" ? mesh.ny*mesh.dy/2.0 : yc  #[0, ny*dy]
    geo.zc =  zc == "unvalued" ? mesh.nz*mesh.dz/2.0 : zc  #[0, nz*dz]
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

"""
    create_box(mesh::Mesh; x1=0, y1=0, z1=0, x2=0, y2=0, z2=0)

Create a cuboid from (x1,y1,z1) to (x2,y2,z2). For example:
```julia
    mesh = FDMesh(nx=10,ny=10,nz=10,dx=1e-9,dy=1e-9,dz=1e-9)
    myBox = create_box(mesh, x1=0, y1=0, z1=0, x2=5e-9, y2=10e-9, z2=10e-9)  
```
Then, one can use functions to set parameters within the box. For example:
```julia
    sim = Sim(mesh)
    set_Ms(sim, myBox, 1e5) 
    add_exch(sim, myBox, 1e-12)
```
Other area are not effected by these functions.
"""
function create_box(mesh::Mesh; x1=0, y1=0, z1=0, x2=0, y2=0, z2=0)
    geo = Box()
    geo.x1,geo.y1,geo.z1 = x1,y1,z1
    geo.x2,geo.y2,geo.z2 = x2,y2,z2
    geo.shape = zeros(Bool, mesh.nxyz)

    for i=1:mesh.nx, j=1:mesh.ny, k=1:mesh.nz
        x = (i-0.5)*mesh.dx
        y = (j-0.5)*mesh.dy
        z = (k-0.5)*mesh.dz
        id = index(i,j,k, mesh.nx, mesh.ny, mesh.nz)
        if x1<=x<=x2 && y1<=y<=y2 && z1<=z<=z2
            geo.shape[id] = true
        end
    end

    return geo
end
