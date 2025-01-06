@inline safe_div(a, b) = b == 0 ? 0.0 : a / b


"""
    The kernel for the Zeeman interaction, works for both the micromagnetic and atomistic model,
where factor = cell_size for the micromagnetic model and factor = 1 for atomistic model.
"""
@kernel function zeeman_kernel!(@Const(m), @Const(h), energy, @Const(mu0_Ms),
                                factor::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds mh = m[j + 1] * h[j + 1] + m[j + 2] * h[j + 2] + m[j + 3] * h[j + 3]
    @inbounds energy[id] = -factor * mu0_Ms[id] * mh
end

"""
Similar to the zeeman_kernel! that this kernel works for both the micromagnetic and atomistic model,
and factor = cell_size for the micromagnetic model and factor = 1 for atomistic model.
"""
@kernel function time_zeeman_kernel!(@Const(m), h, @Const(h_static), energy, @Const(mu0_Ms),
                                     factor::T, fx::T, fy::T,
                                     fz::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * (I - 1)
    @inbounds h[j + 1] = h_static[j + 1] * fx
    @inbounds h[j + 2] = h_static[j + 2] * fy
    @inbounds h[j + 3] = h_static[j + 3] * fz
    @inbounds mh::T = m[j + 1] * h[j + 1] + m[j + 2] * h[j + 2] + m[j + 3] * h[j + 3]
    @inbounds energy[I] = -factor * mu0_Ms[I] * mh
end

"""
The kernel anisotropy_kernel! works for both the micromagnetic and atomistic model, and volume = 1 for atomistic model.
"""
@kernel function anisotropy_kernel!(@Const(m), h, energy, @Const(Ku), axis_x::T, axis_y::T,
                                    axis_z::T, @Const(mu0_Ms),
                                    volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = mu0_Ms[id]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        Ms_inv::T = 2.0 / Ms_local
        @inbounds sa = m[j + 1] * axis_x + m[j + 2] * axis_y + m[j + 3] * axis_z
        @inbounds h[j + 1] = Ku[id] * sa * axis_x * Ms_inv
        @inbounds h[j + 2] = Ku[id] * sa * axis_y * Ms_inv
        @inbounds h[j + 3] = Ku[id] * sa * axis_z * Ms_inv
        @inbounds energy[id] = Ku[id] * (1.0 - sa * sa) * volume
    end
end

"""
The kernel spatial_anisotropy_kernel! works for both the micromagnetic and atomistic model, and volume = 1 for atomistic model.
"""
@kernel function spatial_anisotropy_kernel!(@Const(m), h, energy, @Const(Ku), @Const(axes),
                                            @Const(mu0_Ms),
                                            volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = mu0_Ms[id]

    @inbounds ax = axes[1, I]
    @inbounds ay = axes[2, I]
    @inbounds az = axes[3, I]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        Ms_inv::T = 2.0 / Ms_local
        @inbounds sa = m[j + 1] * axis_x + m[j + 2] * axis_y + m[j + 3] * axis_z
        @inbounds h[j + 1] = Ku[id] * sa * ax * Ms_inv
        @inbounds h[j + 2] = Ku[id] * sa * ay * Ms_inv
        @inbounds h[j + 3] = Ku[id] * sa * az * Ms_inv
        @inbounds energy[id] = Ku[id] * (1.0 - sa * sa) * volume
    end
end

"""
The kernel cubic_anisotropy_kernel! works for both the micromagnetic and atomistic model, and volume = 1 for atomistic model.
"""
@kernel function cubic_anisotropy_kernel!(@Const(m), h, energy, @Const(Kc), a1x::T, a1y::T,
                                          a1z::T, a2x::T, a2y::T, a2z::T, a3x::T, a3y::T,
                                          a3z::T, @Const(mu0_Ms),
                                          volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = mu0_Ms[id]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        Ms_inv::T = 4.0 * Kc[id] / Ms_local
        @inbounds mxp = a1x * m[j + 1] + a1y * m[j + 2] + a1z * m[j + 3]
        @inbounds myp = a2x * m[j + 1] + a2y * m[j + 2] + a2z * m[j + 3]
        @inbounds mzp = a3x * m[j + 1] + a3y * m[j + 2] + a3z * m[j + 3]
        mxp3 = mxp * mxp * mxp
        myp3 = myp * myp * myp
        mzp3 = mzp * mzp * mzp
        @inbounds h[j + 1] = Ms_inv * (mxp3 * a1x + myp3 * a2x + mzp3 * a3x)
        @inbounds h[j + 2] = Ms_inv * (mxp3 * a1y + myp3 * a2y + mzp3 * a3y)
        @inbounds h[j + 3] = Ms_inv * (mxp3 * a1z + myp3 * a2z + mzp3 * a3z)
        @inbounds energy[id] = -Kc[id] * (mxp * mxp3 + myp * myp3 + mzp * mzp3) * volume
    end
end

@kernel function exchange_kernel!(@Const(m), h, energy, @Const(mu0_Ms), @Const(A), dx::T,
                                  dy::T, dz::T, @Const(ngbs),
                                  volume::T) where {T<:AbstractFloat}
    I = @index(Global)

    @inbounds Ms_local = mu0_Ms[I]

    ax::T = 2 / (dx * dx)
    ay::T = 2 / (dy * dy)
    az::T = 2 / (dz * dz)
    nabla = (ax, ax, ay, ay, az, az)

    i = 3 * I - 2
    if Ms_local == T(0)
        @inbounds energy[I] = 0
        @inbounds h[i] = 0
        @inbounds h[i + 1] = 0
        @inbounds h[i + 2] = 0
    else
        fx, fy, fz = T(0), T(0), T(0)
        for j in 1:6
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                k = 3 * id - 2
                @inbounds A_eff = safe_div(2 * A[I] * A[id], A[I] + A[id])
                @inbounds fx += A_eff * nabla[j] * (m[k] - m[i])
                @inbounds fy += A_eff * nabla[j] * (m[k + 1] - m[i + 1])
                @inbounds fz += A_eff * nabla[j] * (m[k + 2] - m[i + 2])
            end
        end
        Ms_inv = 1.0 / (Ms_local)
        @inbounds energy[I] = -0.5 * (fx * m[i] + fy * m[i + 1] + fz * m[i + 2]) * volume
        @inbounds h[i] = fx * Ms_inv
        @inbounds h[i + 1] = fy * Ms_inv
        @inbounds h[i + 2] = fz * Ms_inv
    end
end

@kernel function uniform_exchange_kernel!(@Const(m), h, energy, @Const(mu0_Ms), Ax::T,
                                          Ay::T, Az::T, dx::T, dy::T, dz::T, @Const(ngbs),
                                          volume::T) where {T<:AbstractFloat}
    I = @index(Global)

    @inbounds Ms_local = mu0_Ms[I]

    ax::T = 2 * Ax / (dx * dx)
    ay::T = 2 * Ay / (dy * dy)
    az::T = 2 * Az / (dz * dz)
    nabla = (ax, ax, ay, ay, az, az)

    i = 3 * I - 2
    if Ms_local == T(0)
        @inbounds energy[I] = 0
        @inbounds h[i] = 0
        @inbounds h[i + 1] = 0
        @inbounds h[i + 2] = 0
    else
        fx, fy, fz = T(0), T(0), T(0)
        for j in 1:6
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                k = 3 * id - 2
                @inbounds fx += nabla[j] * (m[k] - m[i])
                @inbounds fy += nabla[j] * (m[k + 1] - m[i + 1])
                @inbounds fz += nabla[j] * (m[k + 2] - m[i + 2])
            end
        end
        Ms_inv = 1.0 / Ms_local
        @inbounds energy[I] = -0.5 * (fx * m[i] + fy * m[i + 1] + fz * m[i + 2]) * volume
        @inbounds h[i] = fx * Ms_inv
        @inbounds h[i + 1] = fy * Ms_inv
        @inbounds h[i + 2] = fz * Ms_inv
    end
end

@kernel function bulkdmi_kernel!(@Const(m), h, energy, @Const(mu0_Ms), Dx::T, Dy::T, Dz::T,
                                 dx::T, dy::T, dz::T, @Const(ngbs),
                                 volume::T) where {T<:AbstractFloat}
    I = @index(Global)

    @inbounds Ms_local = mu0_Ms[I]

    Ds = (Dx / dx, Dx / dx, Dy / dy, Dy / dy, Dz / dz, Dz / dz)
    ax = (T(1), T(-1), T(0), T(0), T(0), T(0))
    ay = (T(0), T(0), T(1), T(-1), T(0), T(0))
    az = (T(0), T(0), T(0), T(0), T(1), T(-1))

    i = 3 * I - 2
    if Ms_local == T(0)
        @inbounds energy[I] = 0
        @inbounds h[i] = 0
        @inbounds h[i + 1] = 0
        @inbounds h[i + 2] = 0
    else
        fx, fy, fz = T(0), T(0), T(0)
        for j in 1:6
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                k = 3 * id - 2
                @inbounds fx += Ds[j] *
                                cross_x(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
                @inbounds fy += Ds[j] *
                                cross_y(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
                @inbounds fz += Ds[j] *
                                cross_z(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
            end
        end
        Ms_inv = 1.0 / Ms_local
        @inbounds energy[I] = -0.5 * (fx * m[i] + fy * m[i + 1] + fz * m[i + 2]) * volume
        @inbounds h[i] = fx * Ms_inv
        @inbounds h[i + 1] = fy * Ms_inv
        @inbounds h[i + 2] = fz * Ms_inv
    end
end

@kernel function spatial_bulkdmi_kernel!(@Const(m), h, energy, @Const(mu0_Ms), @Const(Ds),
                                         dx::T, dy::T, dz::T, @Const(ngbs),
                                         volume::T) where {T<:AbstractFloat}
    I = @index(Global)
    @inbounds Ms_local = mu0_Ms[I]

    Dd = (T(1 / dx), T(1 / dx), T(1 / dy), T(1 / dy), T(1 / dz), T(1 / dz))
    ax = (T(1), T(-1), T(0), T(0), T(0), T(0))
    ay = (T(0), T(0), T(1), T(-1), T(0), T(0))
    az = (T(0), T(0), T(0), T(0), T(1), T(-1))

    i = 3 * I - 2
    if Ms_local == T(0)
        @inbounds energy[I] = 0
        @inbounds h[i] = 0
        @inbounds h[i + 1] = 0
        @inbounds h[i + 2] = 0
    else
        fx, fy, fz = T(0), T(0), T(0)
        for j in 1:6
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                k = 3 * id - 2
                @inbounds D = safe_div(2 * Ds[I] * Ds[id], Ds[I] + Ds[id])
                @inbounds fx += D *
                                Dd[j] *
                                cross_x(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
                @inbounds fy += D *
                                Dd[j] *
                                cross_y(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
                @inbounds fz += D *
                                Dd[j] *
                                cross_z(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
            end
        end
        Ms_inv = 1.0 / Ms_local
        @inbounds energy[I] = -0.5 * (fx * m[i] + fy * m[i + 1] + fz * m[i + 2]) * volume
        @inbounds h[i] = fx * Ms_inv
        @inbounds h[i + 1] = fy * Ms_inv
        @inbounds h[i + 2] = fz * Ms_inv
    end
end

@kernel function interfacial_dmi_kernel!(@Const(m), h, energy, @Const(mu0_Ms), @Const(Ds),
                                         dx::T, dy::T, dz::T, @Const(ngbs),
                                         volume::T) where {T<:AbstractFloat}
    I = @index(Global)
    @inbounds Ms_local = mu0_Ms[I]

    Dd = (T(1 / dx), T(1 / dx), T(1 / dy), T(1 / dy))
    ax = (T(0), T(0), T(-1), T(1))
    ay = (T(1), T(-1), T(0), T(0))
    az = (T(0), T(0), T(0), T(0))

    i = 3 * I - 2
    if Ms_local == T(0)
        @inbounds energy[I] = 0
        @inbounds h[i] = 0
        @inbounds h[i + 1] = 0
        @inbounds h[i + 2] = 0
    else
        fx, fy, fz = T(0), T(0), T(0)
        for j in 1:4
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                k = 3 * id - 2
                @inbounds D = safe_div(2 * Ds[I] * Ds[id], Ds[I] + Ds[id])
                @inbounds fx += D *
                                Dd[j] *
                                cross_x(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
                @inbounds fy += D *
                                Dd[j] *
                                cross_y(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
                @inbounds fz += D *
                                Dd[j] *
                                cross_z(ax[j], ay[j], az[j], m[k], m[k + 1], m[k + 2])
            end
        end
        Ms_inv = 1.0 / Ms_local
        @inbounds energy[I] = -0.5 * (fx * m[i] + fy * m[i + 1] + fz * m[i + 2]) * volume
        @inbounds h[i] = fx * Ms_inv
        @inbounds h[i + 1] = fy * Ms_inv
        @inbounds h[i + 2] = fz * Ms_inv
    end
end

@kernel function stochastic_field_kernel!(@Const(m), h, energy, @Const(mu0_Ms), @Const(eta),
                                          @Const(Temp), factor::T,
                                          volume::T) where {T<:AbstractFloat}
    I = @index(Global)

    j = 3 * (I - 1)
    @inbounds Ms_local = mu0_Ms[I]
    @inbounds T_local = Temp[I]

    if Ms_local > 0
        @inbounds scale = sqrt(factor * T_local / Ms_local)
        @inbounds h[j + 1] = eta[j + 1] * scale
        @inbounds h[j + 2] = eta[j + 2] * scale
        @inbounds h[j + 3] = eta[j + 3] * scale
        @inbounds energy[I] = -Ms_local *
                              volume *
                              (m[j + 1] * h[j + 1] +
                               m[j + 2] * h[j + 2] +
                               m[j + 3] * h[j + 3])
    else
        @inbounds energy[I] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    end
end

@kernel function interlayer_exch_kernel!(@Const(m), h, energy, @Const(mu0_Ms), J::T,
                                         K1::Int32, K2::Int32, nx::Int32, ny::Int32, dz::T,
                                         volume::T) where {T<:AbstractFloat}
    i, j = @index(Global, NTuple)

    id1 = (K1 - 1) * nx * ny + (j - 1) * nx + i
    id2 = (K2 - 1) * nx * ny + (j - 1) * nx + i

    k1 = 3 * id1 - 2
    k2 = 3 * id2 - 2
    @inbounds mbx = m[k1]
    @inbounds mby = m[k1 + 1]
    @inbounds mbz = m[k1 + 2]

    @inbounds mtx = m[k2]
    @inbounds mty = m[k2 + 1]
    @inbounds mtz = m[k2 + 2]

    @inbounds Ms1 = mu0_Ms[id1]
    @inbounds Ms2 = mu0_Ms[id2]
    if Ms1 > 0 && Ms2 > 0
        Ms_inv = J / (Ms1 * dz)
        @inbounds h[k1] = Ms_inv * mtx
        @inbounds h[k1 + 1] = Ms_inv * mty
        @inbounds h[k1 + 2] = Ms_inv * mtz
        @inbounds energy[id1] = -0.5 *
                                (h[k1] * mbx + h[k1 + 1] * mby + h[k1 + 2] * mbz) *
                                volume

        Ms_inv = J / (Ms2 * dz)
        @inbounds h[k2] = Ms_inv * mbx
        @inbounds h[k2 + 1] = Ms_inv * mby
        @inbounds h[k2 + 2] = Ms_inv * mbz
        @inbounds energy[id2] = -0.5 *
                                (h[k2] * mtx + h[k2 + 1] * mty + h[k2 + 2] * mtz) *
                                volume
    end
end

@kernel function interlayer_dmi_kernel!(@Const(m), h, energy, @Const(mu0_Ms), Dx::T, Dy::T,
                                        Dz::T, K1::Int32, K2::Int32, nx::Int32, ny::Int32,
                                        dz::T, volume::T) where {T<:AbstractFloat}
    i, j = @index(Global, NTuple)

    id1 = (K1 - 1) * nx * ny + (j - 1) * nx + i
    id2 = (K2 - 1) * nx * ny + (j - 1) * nx + i

    k1 = 3 * id1 - 2
    k2 = 3 * id2 - 2
    @inbounds mbx = m[k1]
    @inbounds mby = m[k1 + 1]
    @inbounds mbz = m[k1 + 2]

    @inbounds mtx = m[k2]
    @inbounds mty = m[k2 + 1]
    @inbounds mtz = m[k2 + 2]

    @inbounds Ms1 = mu0_Ms[id1]
    @inbounds Ms2 = mu0_Ms[id2]
    if Ms1 > 0 && Ms2 > 0
        Ms_inv = 1.0 / (Ms1 * dz)
        @inbounds h[k1] = Ms_inv * cross_x(Dx, Dy, Dz, mtx, mty, mtz)
        @inbounds h[k1 + 1] = Ms_inv * cross_y(Dx, Dy, Dz, mtx, mty, mtz)
        @inbounds h[k1 + 2] = Ms_inv * cross_z(Dx, Dy, Dz, mtx, mty, mtz)
        @inbounds energy[id1] = -0.5 *
                                (h[k1] * mbx + h[k1 + 1] * mby + h[k1 + 2] * mbz) *
                                volume

        Ms_inv = -1.0 / (Ms2 * dz)
        @inbounds h[k2] = Ms_inv * cross_x(Dx, Dy, Dz, mbx, mby, mbz)
        @inbounds h[k2 + 1] = Ms_inv * cross_y(Dx, Dy, Dz, mbx, mby, mbz)
        @inbounds h[k2 + 2] = Ms_inv * cross_z(Dx, Dy, Dz, mbx, mby, mbz)
        @inbounds energy[id2] = -0.5 *
                                (h[k2] * mtx + h[k2 + 1] * mty + h[k2 + 2] * mtz) *
                                volume
    end
end
