using LinearAlgebra

@kernel function cartesian2spherical_kernel!(@Const(m), spherical)
    id = @index(Global)
    j = 3 * (id - 1)
    norm = sqrt(m[j+1]*m[j+1] + m[j+2]*m[j+2] + m[j+3]*m[j+3])
    @inbounds spherical[j+1] = norm
    @inbounds spherical[j+2] = acos(m[j+3]/(norm+1e-30))
    @inbounds spherical[j+3] = atan(m[j+2],m[j+1])
end

@kernel function spherical2cartesian_kernel!(@Const(spherical), cartesian)
    id = @index(Global)
    j = 3 * (id - 1)
    r, theta, psi = spherical[j+1], spherical[j+2], spherical[j+3]
    @inbounds cartesian[j+1] = r * sin(theta) * cos(psi)
    @inbounds cartesian[j+2] = r * sin(theta) * sin(psi)
    @inbounds cartesian[j+3] = r * cos(theta)
end

@kernel function inner_product_kernel!(@Const(x1), @Const(x2), res)
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds res[id] = x1[j+1]*x2[j+1] + x1[j+2]*x2[j+2] + x1[j+3]*x2[j+3]
end

@kernel function slerp_kernel!(@Const(x1), @Const(x2), t, res)
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds inner = x1[j+1]*x2[j+1] + x1[j+2]*x2[j+2] + x1[j+3]*x2[j+3]

    theta = acos(inner)
    if !(theta == 0.)
        sin_theta = sin(theta)
        a = sin((1-t)*theta) / sin_theta
        b = sin(t*theta) / sin_theta
        @inbounds res[j+1] = a * x1[j+1] + b * x2[j+1]
        @inbounds res[j+2] = a * x1[j+2] + b * x2[j+2]
        @inbounds res[j+3] = a * x1[j+3] + b * x2[j+3]
    end
end

"""
Bessarab P F, et al. Computer Physics Communications, 2015, 196: 335-347.
Appendix B.
"""
@kernel function compute_distance_kernel!(@Const(x1), @Const(x2), distance)
    id = @index(Global)
    j = 3 * (id - 1)
    dot_p = dot_product(x1[j+1], x1[j+2], x1[j+3], x2[j+1], x2[j+2], x2[j+3])
    px, py, pz = cross_product(x1[j+1], x1[j+2], x1[j+3], x2[j+1], x2[j+2], x2[j+3])
    n_cross = hypot(px,py,pz)
    distance[id] = atan(n_cross, dot_p)
end

@kernel function reduce_tangent_kernel!(field, @Const(tangents))
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds f = dot_product(field[j+1],field[j+2],field[j+3],tangents[j+1],tangents[j+2],tangents[j+3])
    @inbounds field[j+1] -= f*tangents[j+1]
    @inbounds field[j+2] -= f*tangents[j+2]
    @inbounds field[j+3] -= f*tangents[j+3]
end