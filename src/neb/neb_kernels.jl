using LinearAlgebra

@kernel function cartesian2spherical_kernel!(@Const(m), spherical)
    id = @index(Global)
    j = 3 * (id - 1)
    norm = sqrt(m[j + 1] * m[j + 1] + m[j + 2] * m[j + 2] + m[j + 3] * m[j + 3])
    @inbounds spherical[j + 1] = norm
    @inbounds spherical[j + 2] = acos(m[j + 3] / (norm + 1e-30))
    @inbounds spherical[j + 3] = atan(m[j + 2], m[j + 1])
end

@kernel function spherical2cartesian_kernel!(@Const(spherical), cartesian)
    id = @index(Global)
    j = 3 * (id - 1)
    r, theta, psi = spherical[j + 1], spherical[j + 2], spherical[j + 3]
    @inbounds cartesian[j + 1] = r * sin(theta) * cos(psi)
    @inbounds cartesian[j + 2] = r * sin(theta) * sin(psi)
    @inbounds cartesian[j + 3] = r * cos(theta)
end

@kernel function inner_product_kernel!(@Const(x1), @Const(x2), res)
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds res[id] = x1[j + 1] * x2[j + 1] +
                        x1[j + 2] * x2[j + 2] +
                        x1[j + 3] * x2[j + 3]
end

@kernel function slerp_kernel!(@Const(x1), @Const(x2), t, res)
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds inner = x1[j + 1] * x2[j + 1] + x1[j + 2] * x2[j + 2] + x1[j + 3] * x2[j + 3]

    theta = acos(inner)
    if !(theta == 0.0)
        sin_theta = sin(theta)
        a = sin((1 - t) * theta) / sin_theta
        b = sin(t * theta) / sin_theta
        @inbounds res[j + 1] = a * x1[j + 1] + b * x2[j + 1]
        @inbounds res[j + 2] = a * x1[j + 2] + b * x2[j + 2]
        @inbounds res[j + 3] = a * x1[j + 3] + b * x2[j + 3]
    end
end

"""
Bessarab P F, et al. Computer Physics Communications, 2015, 196: 335-347.
Appendix B.
"""
@kernel function compute_distance_kernel!(@Const(m1), @Const(m2), distance)
    id = @index(Global)
    j = 3 * id - 2
    @inbounds mx1 = m1[j]
    @inbounds my1 = m1[j + 1]
    @inbounds mz1 = m1[j + 2]
    @inbounds mx2 = m2[j]
    @inbounds my2 = m2[j + 1]
    @inbounds mz2 = m2[j + 2]
    a = -my2 * mz1 + my1 * mz2  #m1xm2, x
    b = mx2 * mz1 - mx1 * mz2   #m1xm2, y
    c = -mx2 * my1 + mx1 * my2  #m1xm2, z
    mm = mx1 * mx2 + my1 * my2 + mz1 * mz2
    d = sqrt(a * a + b * b + c * c)
    @inbounds distance[id] = atan(d, mm)
end

#compute h = h - (m.h)m
@kernel function reduce_tangent_kernel!(field, @Const(m))
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds f = dot_product(field[j + 1], field[j + 2], field[j + 3], m[j + 1], m[j + 2],
                              m[j + 3])
    @inbounds field[j + 1] -= f * m[j + 1]
    @inbounds field[j + 2] -= f * m[j + 2]
    @inbounds field[j + 3] -= f * m[j + 3]
end

@kernel function compute_tangents_kernel!(t, m, left_m, right_m, energy, N::Int64,
                                          dof::Int64)
    id = @index(Global)

    for n in 1:N
        k = (n - 1) * dof + id
        E1 = energy[n]
        E2 = energy[n + 1]
        E3 = energy[n + 2]
        dEmax = max(abs(E3 - E2), abs(E2 - E1))
        dEmin = min(abs(E3 - E2), abs(E2 - E1))

        m_kp = n == N ? right_m[id] : m[k + dof]
        m_k = m[k]
        m_km = n == 1 ? left_m[id] : m[k - dof]
        tip = m_kp - m_k
        tim = m_k - m_km

        if (E1 > E2) && (E2 > E3)
            t[k] = tim
        elseif (E3 > E2) && (E2 > E1)
            t[k] = tip
        elseif E3 > E1
            t[k] = dEmax * tip + dEmin * tim
        elseif E3 < E1
            t[k] = dEmin * tip + dEmax * tim
        else
            t[k] = tim + tip
        end
    end
end
