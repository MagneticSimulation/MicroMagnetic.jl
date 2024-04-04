@kernel function zeeman_kernel!(@Const(m), @Const(h), energy, @Const(Ms), volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds energy[id] = -mu_0 * Ms[id] * volume * (m[j+1] * h[j+1] + m[j+2] * h[j+2] + m[j+3] * h[j+3])
end


@kernel function time_zeeman_kernel!(@Const(m), h, @Const(h_static), energy, @Const(Ms), volume::T, fx::T, fy::T, fz::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds h[j+1] = h_static[j+1] * fx
    @inbounds h[j+2] = h_static[j+2] * fy
    @inbounds h[j+3] = h_static[j+3] * fz
    @inbounds energy[id] = -mu_0 * Ms[id] * volume * (m[j+1] * h[j+1] + m[j+2] * h[j+2] + m[j+3] * h[j+3])
end


@kernel function anisotropy_kernel!(@Const(m), h, energy, @Const(Ku), axis_x::T, axis_y::T, axis_z::T, @Const(Ms), volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = Ms[id]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j+1] = 0
        @inbounds h[j+2] = 0
        @inbounds h[j+3] = 0
    else
        Ms_inv::T = 2.0 / (mu_0 * Ms_local)
        @inbounds sa = m[j+1] * axis_x + m[j+2] * axis_y + m[j+3] * axis_z
        @inbounds h[j+1] = Ku[id] * sa * axis_x * Ms_inv
        @inbounds h[j+2] = Ku[id] * sa * axis_y * Ms_inv
        @inbounds h[j+3] = Ku[id] * sa * axis_z * Ms_inv
        @inbounds energy[id] = Ku[id] * (1.0 - sa * sa) * volume
    end

end

@kernel function cubic_anisotropy_kernel!(@Const(m), h, energy, @Const(Kc), a1x::T, a1y::T, a1z::T,
    a2x::T, a2y::T, a2z::T, a3x::T, a3y::T, a3z::T, @Const(Ms), volume::T) where {T<:AbstractFloat}

    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = Ms[id]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j+1] = 0
        @inbounds h[j+2] = 0
        @inbounds h[j+3] = 0
    else
        Ms_inv::T = 4.0 * Kc[id] / (mu_0 * Ms_local)
        @inbounds mxp::T = a1x * m[j+1] + a1y * m[j+2] + a1z * m[j+3]
        @inbounds myp::T = a2x * m[j+1] + a2y * m[j+2] + a2z * m[j+3]
        @inbounds mzp::T = a3x * m[j+1] + a3y * m[j+2] + a3z * m[j+3]
        mxp3::T = mxp * mxp * mxp
        myp3::T = myp * myp * myp
        mzp3::T = mzp * mzp * mzp
        @inbounds h[j+1] = Ms_inv * (mxp3 * a1x + myp3 * a2x + mzp3 * a3x)
        @inbounds h[j+2] = Ms_inv * (mxp3 * a1y + myp3 * a2y + mzp3 * a3y)
        @inbounds h[j+3] = Ms_inv * (mxp3 * a1z + myp3 * a2z + mzp3 * a3z)
        @inbounds energy[id] = -Kc[id] * (mxp * mxp3 + myp * myp3 + mzp * mzp3) * volume
    end
end


@kernel function exchange_kernel!(@Const(m), h, energy, @Const(Ms), @Const(A),
    dx::T, dy::T, dz::T, @Const(ngbs), volume::T) where {T<:AbstractFloat}

    I = @index(Global)

    @inbounds Ms_local = Ms[I]

    ax::T = 2 / (dx * dx)
    ay::T = 2 / (dy * dy)
    az::T = 2 / (dz * dz)
    nabla = (ax, ax, ay, ay, az, az)

    i = 3 * I - 2
    if Ms_local == T(0)
        @inbounds energy[I] = 0
        @inbounds h[i] = 0
        @inbounds h[i+1] = 0
        @inbounds h[i+2] = 0
    else
        fx, fy, fz = T(0), T(0), T(0)
        for j = 1:6
            @inbounds id = ngbs[j, I]
            if id > 0 && Ms_local > 0
                k = 3 * id - 2
                @inbounds fx += A[I] * nabla[j] * (m[k] - m[i])
                @inbounds fy += A[I] * nabla[j] * (m[k+1] - m[i+1])
                @inbounds fz += A[I] * nabla[j] * (m[k+2] - m[i+2])
            end
        end
        Ms_inv = 1.0 / (Ms_local * mu_0)
        @inbounds energy[I] = -0.5 * (fx * m[i] + fy * m[i+1] + fz * m[i+2]) * volume
        @inbounds h[i] = fx * Ms_inv
        @inbounds h[i+1] = fy * Ms_inv
        @inbounds h[i+2] = fz * Ms_inv
    end
end