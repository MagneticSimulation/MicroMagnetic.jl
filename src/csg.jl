
# constructive solid geometry
# We implement a simple version that cotains Plane, Sphere, Cylinder, Box and Torus.

using LinearAlgebra
import Base: +, -, *

export +, -, *, Plane, Sphere, Cylinder, Box, Torus
export CSGShape

abstract type CSGShape end
const Shape = CSGShape

struct UnionShape <: CSGShape
    left::CSGShape
    right::CSGShape
end

struct IntersectionShape <: CSGShape
    left::CSGShape
    right::CSGShape
end

struct DifferenceShape <: CSGShape
    left::CSGShape
    right::CSGShape
end

function normalize(x::Tuple{Real,Real,Real})
    length = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    return x ./ length
end

"""
    Plane(; point::Tuple=(0.0, 0.0, 0.0), normal::Tuple=(0.0, 0.0, 1.0))

Create a Plane shape.

Parameters
----------
- point: A point on the plane
- normal: Normal vector of the plane (will be normalized automatically)

Examples
--------
```julia
# Create a plane at origin with normal along z-axis
plane = Plane()

# Create a plane at (1, 2, 3) with normal (1, 1, 1)
plane = Plane(point=(1, 2, 3), normal=(1, 1, 1))
```
"""
struct Plane <: CSGShape
    point::Tuple{Real,Real,Real} # A point in the plane
    normal::Tuple{Real,Real,Real} # the normal of the plane

    function Plane(; point::Tuple=(0.0, 0.0, 0.0), normal::Tuple=(0.0, 0.0, 1.0))
        return new(point, normalize(normal))
    end
end

"""
    Sphere(; center::Tuple=(0.0, 0.0, 0.0), radius=1.0)

Create a Sphere shape.

Parameters
----------
- center: Center of the sphere
- radius: Radius of the sphere

Examples
--------
```julia
# Create a sphere at origin with radius 1.0
sphere = Sphere()

# Create a sphere at (2, 0, 0) with radius 5.0
sphere = Sphere(center=(2, 0, 0), radius=5.0)
```
"""
struct Sphere <: CSGShape # A sphere
    center::Tuple{Real,Real,Real}
    radius::Float64

    function Sphere(; center::Tuple=(0.0, 0.0, 0.0), radius=1.0)
        return new(center, radius)
    end
end

"""
    Cylinder(; center=(0.0, 0.0, 0.0), radius=1.0, height=Inf, normal=(0.0, 0.0, 1.0))

Create a Cylinder shape.

Parameters
----------
- center: Center of the cylinder
- radius: Radius of the cylinder
- height: Height of the cylinder (default: Inf, meaning infinite length)
- normal: Direction of the cylinder axis (will be normalized automatically)

Examples
--------
```julia
# Create an infinite cylinder along z-axis
cylinder = Cylinder()

# Create a finite cylinder along y-axis
cylinder = Cylinder(center=(0, 0, 0), radius=2.0, height=10.0, normal=(0, 1, 0))
```
"""
struct Cylinder <: Shape
    center::Tuple{Real,Real,Real}
    radius::Float64
    height::Float64
    normal::Tuple{Real,Real,Real} #should be normalized

    function Cylinder(; center::Tuple=(0.0, 0.0, 0.0), radius=1.0, height=Inf,
                      normal::Tuple=(0.0, 0.0, 1.0))
        n = normalize(normal)
        return new(center, radius, height, n)
    end
end

"""
    Box(; center::Tuple=(0.0, 0.0, 0.0), size::Tuple=(1.0, 1.0, 1.0), theta=0.0)

Create a Box shape.

Parameters
----------
- center: Center of the box
- size: Size of the box in (x, y, z) directions
- theta: Rotation angle around z-axis (in radians)

Examples
--------
```julia
# Create a unit box at origin
box = Box()

# Create a larger box with rotation
box = Box(center=(1, 1, 1), size=(2, 3, 4), theta=pi/4)
```
"""
struct Box <: CSGShape
    center::Tuple{Real,Real,Real}
    size::Tuple{Real,Real,Real}
    theta::Float64 # the rotated angle
    function Box(; center::Tuple=(0.0, 0.0, 0.0), size::Tuple=(1.0, 1.0, 1.0), theta=0.0)
        return new(center, size, theta)
    end
end

"""
    Torus(; center::Tuple=(0.0, 0.0, 0.0), R=1.0, r=0.2)

Create a Torus shape.

Parameters
----------
- center: Center of the torus
- R: Major radius (distance from center of torus to center of tube)
- r: Minor radius (radius of the tube)

Examples
--------
```julia
# Create a torus with major radius 1.0 and minor radius 0.2
 torus = Torus()

# Create a torus with custom dimensions
 torus = Torus(center=(0, 0, 0), R=5.0, r=1.0)
```
"""
struct Torus <: CSGShape
    center::Tuple{Real,Real,Real} # Center of the torus
    R::Float64 # Major radius
    r::Float64 # Minor radius

    # Constructor for Torus struct
    function Torus(; center::Tuple=(0.0, 0.0, 0.0), R=1.0, r=0.2)
        return new(center, R, r)
    end
end

+(a::CSGShape, b::CSGShape) = UnionShape(a, b)
*(a::CSGShape, b::CSGShape) = IntersectionShape(a, b)
-(a::CSGShape, b::CSGShape) = DifferenceShape(a, b)

# Check if the given point is inside the halfspace
function point_inside_halfspace(point::Tuple{Real,Real,Real}, plane::Plane)
    dot_product = dot(point .- plane.point, plane.normal)
    return dot_product >= 0
end

# Check if the given point is inside the sphere
function point_inside_sphere(point::Tuple{Real,Real,Real}, sphere::Sphere)
    distance = sqrt(sum((point[i] - sphere.center[i])^2 for i in 1:3))
    return distance <= sphere.radius
end

# Check if the given point is inside the Cylinder
function point_inside_cylinder(point::Tuple{Real,Real,Real}, cylinder::Cylinder)
    d = cross_product(point .- cylinder.center, cylinder.normal)
    if dot(d, d) > cylinder.radius^2
        return false
    end

    # Check the z direction
    height_difference = abs(point[3] - cylinder.center[3])
    if height_difference > cylinder.height / 2.0
        return false
    end

    return true
end

# Check if the given point is inside the Box
function point_inside_box(point::Tuple{Real,Real,Real}, box::Box)
    local_point = (point[1] - box.center[1], point[2] - box.center[2],
                   point[3] - box.center[3])

    # Rotate the point back to the original coordinate system
    rotated_x = local_point[1] * cos(-box.theta) - local_point[2] * sin(-box.theta)
    rotated_y = local_point[1] * sin(-box.theta) + local_point[2] * cos(-box.theta)
    rotated_point = (rotated_x, rotated_y, local_point[3])

    # Check if the rotated point is inside the Box
    half_width = box.size[1] / 2
    half_height = box.size[2] / 2
    half_depth = box.size[3] / 2

    if -half_width <= rotated_point[1] <= half_width &&
       -half_height <= rotated_point[2] <= half_height &&
       -half_depth <= rotated_point[3] <= half_depth
        return true
    else
        return false
    end
end

# see https://math.stackexchange.com/questions/4380905/how-can-i-determine-if-a-point-x-y-z-is-within-a-torus-r-r
function point_inside_torus(point::Tuple{Real,Real,Real}, torus::Torus)
    p = point .- torus.center

    return (torus.R - sqrt(p[1]^2 + p[2]^2))^2 + p[3]^2 <= torus.r^2
end


function inside(shape::CSGShape, point::Tuple{Real,Real,Real})
    if shape isa Plane
        return point_inside_halfspace(point, shape)
    elseif shape isa Sphere
        return point_inside_sphere(point, shape)
    elseif shape isa Cylinder
        return point_inside_cylinder(point, shape)
    elseif shape isa Box
        return point_inside_box(point, shape)
    elseif shape isa Torus
        return point_inside_torus(point, shape)
    elseif shape isa UnionShape
        return inside(shape.left, point) || inside(shape.right, point)
    elseif shape isa IntersectionShape
        return inside(shape.left, point) && inside(shape.right, point)
    elseif shape isa DifferenceShape
        return inside(shape.left, point) && !inside(shape.right, point)
    else
        error("Unsupported CSGShape type: $(typeof(shape))")
    end
end

"""
    save_vtk(mesh::Mesh, shape::CSGShape, fname::String)

Save the shape to vtk. 

```julia
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=50)
    t1 = Torus(R = 60e-9, r=20e-9)
    save_vtk(mesh, t1, "torus")
```
"""
function save_vtk(mesh::Mesh, shape::CSGShape, fname::String)
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    xyz = zeros(Float32, 3, nx + 1, ny + 1, nz + 1)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz

    scale_factor = 10^floor(log10(dx))
    for k in 1:(nz + 1), j in 1:(ny + 1), i in 1:(nx + 1)
        xyz[1, i, j, k] = (i - 1 - nx / 2) * dx / scale_factor
        xyz[2, i, j, k] = (j - 1 - ny / 2) * dy / scale_factor
        xyz[3, i, j, k] = (k - 1 - nz / 2) * dz / scale_factor
    end
    vtk = vtk_grid(fname, xyz)

    data = zeros(Int32, nx, ny, nz)

    for k in 1:nz, j in 1:ny, i in 1:nx
        x = (i - 0.5 - nx / 2) * dx
        y = (j - 0.5 - ny / 2) * dy
        z = (k - 0.5 - nz / 2) * dz
        if inside(shape, (x, y, z))
            data[i, j, k] = 1
        end
    end

    vtk_cell_data(vtk, data, "shape")

    if scale_factor != 1.0
        vtk["scale_factor", VTKFieldData()] = string(scale_factor)
    end

    return vtk_save(vtk)
end
