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

