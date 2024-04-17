
# constructive solid geometry
# We implement a simple version that cotains Plane, Sphere, Cylinder, Box and Torus.

using LinearAlgebra
import Base: +, -, *

export +, -, *, Plane, Sphere, Cylinder, Box, Torus

abstract type Shape end
const UnionOp = :Union
const IntersectionOp = :Intersection
const DifferenceOp = :Difference

function normalize(x::Tuple{Real,Real,Real})
    length = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    return x ./ length
end

struct Plane <: Shape
    point::Tuple{Real,Real,Real} # A point in the plane
    normal::Tuple{Real,Real,Real} # the normal of the plane

    function Plane(; point::Tuple=(0.0, 0.0, 0.0), normal::Tuple=(0.0, 0.0, 1.0))
        return new(point, normalize(normal))
    end
end

struct Sphere <: Shape # A sphere
    center::Tuple{Real,Real,Real}
    radius::Float64

    function Sphere(; center::Tuple=(0.0, 0.0, 0.0), radius=1.0)
        return new(center, radius)
    end
end

"""
    Cylinder(;center = (0.0, 0.0, 0.0), radius = 1.0, height= Inf, normal = (0.0, 0.0, 1.0))

Create a Cylinder shape. 
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
    Box(;center::Tuple=(0.0, 0.0, 0.0), sides::Tuple=(1.0,1.0,1.0), theta=0.0)

Create a Box shape.
"""
struct Box <: Shape
    center::Tuple{Real,Real,Real}
    sides::Tuple{Real,Real,Real}
    theta::Float64 # the rotated angle
    function Box(; center::Tuple=(0.0, 0.0, 0.0), sides::Tuple=(1.0, 1.0, 1.0), theta=0.0)
        return new(center, sides, theta)
    end
end

"""
    Torus(; center::Tuple=(0.0, 0.0, 0.0), R=1.0, r=0.2)

Create a Torus shape.
"""
struct Torus <: Shape
    center::Tuple{Real,Real,Real} # Center of the torus
    R::Float64 # Major radius
    r::Float64 # Minor radius

    # Constructor for Torus struct
    function Torus(; center::Tuple=(0.0, 0.0, 0.0), R=1.0, r=0.2)
        return new(center, R, r)
    end
end

struct CSGNode
    operation::Symbol
    left::Union{CSGNode,Shape}
    right::Union{CSGNode,Shape}
end

# compute the union of two shapes
function +(shape1::Union{CSGNode,Shape}, shape2::Union{CSGNode,Shape})
    return CSGNode(UnionOp, shape1, shape2)
end

# compute the intersection of two shapes
function *(shape1::Union{CSGNode,Shape}, shape2::Union{CSGNode,Shape})
    return CSGNode(IntersectionOp, shape1, shape2)
end

# compute the difference of two shapes, return shape1 - shape2
function -(shape1::Union{CSGNode,Shape}, shape2::Union{CSGNode,Shape})
    return CSGNode(DifferenceOp, shape1, shape2)
end

# Check if the given point is inside the halfspace
function point_in_halfspace(point::Tuple{Real,Real,Real}, plane::Plane)
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
    half_width = box.sides[1] / 2
    half_height = box.sides[2] / 2
    half_depth = box.sides[3] / 2

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

function point_inside_shape(point::Tuple{Real,Real,Real}, node::CSGNode)
    left_inside = point_inside_shape(point, node.left)

    right_inside = point_inside_shape(point, node.right)

    if node.operation == UnionOp
        return left_inside || right_inside
    elseif node.operation == IntersectionOp
        return left_inside && right_inside
    elseif node.operation == DifferenceOp
        return left_inside && !right_inside
    else
        error("Unsupported boolean operation")
    end
end

function point_inside_shape(point::Tuple{Real,Real,Real}, shape::Shape)
    if isa(shape, Plane)
        return point_in_halfspace(point, shape)
    elseif isa(shape, Sphere)
        return point_inside_sphere(point, shape)
    elseif isa(shape, Cylinder)
        return point_inside_cylinder(point, shape)
    elseif isa(shape, Box)
        return point_inside_box(point, shape)
    elseif isa(shape, Torus)
        return point_inside_torus(point, shape)
    else
        error("Unsupported shape type")
    end
end

"""
    save_vtk(mesh::Mesh, shape::Union{CSGNode,Shape}, fname::String)

Save the shape to vtk. 

```julia
    mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=50)
    t1 = Torus(R = 60e-9, r=20e-9)
    save_vtk(mesh, t1, "torus")
```
"""
function save_vtk(mesh::Mesh, shape::Union{CSGNode,Shape}, fname::String)
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    xyz = zeros(Float32, 3, nx + 1, ny + 1, nz + 1)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    for k in 1:(nz + 1), j in 1:(ny + 1), i in 1:(nx + 1)
        xyz[1, i, j, k] = (i - 1 - nx / 2) * dx
        xyz[2, i, j, k] = (j - 1 - ny / 2) * dy
        xyz[3, i, j, k] = (k - 1 - nz / 2) * dz
    end
    vtk = vtk_grid(fname, xyz)

    data = zeros(Int32, nx, ny, nz)

    for k in 1:nz, j in 1:ny, i in 1:nx
        x = (i - 0.5 - nx / 2) * dx
        y = (j - 0.5 - ny / 2) * dy
        z = (k - 0.5 - nz / 2) * dz
        if point_inside_shape((x, y, z), shape)
            data[i, j, k] = 1
        end
    end

    vtk_cell_data(vtk, data, "shape")

    return vtk_save(vtk)
end
