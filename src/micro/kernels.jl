@inline safe_div(a, b) = b == 0 ? 0.0 : a / b

@kernel function spatiotemporal_kernel!(output, dx::T, dy::T, dz::T, x0::T, y0::T, z0::T,
                                        t::T, f::Function) where {T<:AbstractFloat}
    i, j, k = @index(Global, NTuple)
    x::T = x0 + (i-0.5)*dx
    y::T = y0 + (j-0.5)*dy
    z::T = z0 + (k-0.5)*dz
    @inbounds output[i, j, k] = f(x, y, z, t)
end

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
    The kernel for the Zeeman interaction, works for both the FE micromagnetic model.
"""
@kernel function zeeman_fe_kernel!(@Const(m), @Const(h), energy, @Const(L_mu),
                                   factor::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds mh = m[j + 1] * h[j + 1] * L_mu[j+1] +
                   m[j + 2] * h[j + 2] * L_mu[j+2] +
                   m[j + 3] * h[j + 3] * L_mu[j+3]
    @inbounds energy[id] = -factor * mh
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

    @inbounds ax = axes[j+1]
    @inbounds ay = axes[j+2]
    @inbounds az = axes[j+3]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        Ms_inv::T = 2.0 / Ms_local
        @inbounds sa = m[j + 1] * ax + m[j + 2] * ay + m[j + 3] * az
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

"""
The kernel hexagonal_anisotropy_kernel! works for both the micromagnetic and atomistic model, and volume = 1 for atomistic model.
"""
@kernel function hexagonal_anisotropy_kernel!(@Const(m), h, energy, K1::T, K2::T, K3::T,
                                              @Const(mu0_Ms),
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
        Ms_inv::T = 1.0 / Ms_local
        @inbounds mx = m[j + 1]
        @inbounds my = m[j + 2]
        @inbounds mz = m[j + 3]
        @inbounds h[j + 1] = -6*K3*Ms_inv*(mx^5-10*mx^3*my^2+5*mx*my^4)
        @inbounds h[j + 2] = -6*K3*Ms_inv*(-5*mx^4*my+10*mx^2*my^3-my^5)
        @inbounds h[j + 3] = 2*mz*Ms_inv*(K1 + 2*K2*(1-mz*mz))
        @inbounds energy[id] = (K1*(1-mz*mz) +
                                K2*(1-mz*mz)^2 +
                                K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)) * volume
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
                                          @Const(Temp), base_T::T, factor::T,
                                          volume::T) where {T<:AbstractFloat}
    I = @index(Global)

    j = 3 * (I - 1)
    @inbounds Ms_local = mu0_Ms[I]
    @inbounds T_local = Temp[I] + base_T

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

@kernel function interlayer_exch_kernel!(@Const(m), h, energy, @Const(mu0_Ms), @Const(Js),
                                         K1::Int32, K2::Int32, nx::Int32, ny::Int32, dz::T,
                                         volume::T) where {T<:AbstractFloat}
    i, j = @index(Global, NTuple)

    id1 = (K1 - 1) * nx * ny + (j - 1) * nx + i
    id2 = (K2 - 1) * nx * ny + (j - 1) * nx + i
    id = (j - 1) * nx + i

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
    @inbounds J = Js[id]
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

"""
The kernel sahe_torque_kernel! compute the effective field defined as 
        (1/gamma)*(beta*sigma + m x sigma)
and sigma = c1 - (m.c2)^2 c3
"""
@kernel function sahe_torque_kernel!(@Const(m), h, @Const(mu0_Ms), gamma::T, beta::T,
                                     @Const(c1), c2x::T, c2y::T, c2z::T,
                                     @Const(c3)) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = mu0_Ms[id]

    if Ms_local == 0.0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        @inbounds sa = m[j + 1] * c2x + m[j + 2] * c2y + m[j + 3] * c2z
        @inbounds sx = (c1[j + 1] - sa * sa * c3[j + 1]) / gamma
        @inbounds sy = (c1[j + 2] - sa * sa * c3[j + 2]) / gamma
        @inbounds sz = (c1[j + 3] - sa * sa * c3[j + 3]) / gamma
        @inbounds h[j + 1] = beta*sx + cross_x(m[j + 1], m[j + 2], m[j + 3], sx, sy, sz)
        @inbounds h[j + 2] = beta*sy + cross_y(m[j + 1], m[j + 2], m[j + 3], sx, sy, sz)
        @inbounds h[j + 3] = beta*sz + cross_z(m[j + 1], m[j + 2], m[j + 3], sx, sy, sz)
    end
end

"""
The kernel df_torque_kernel! compute the effective field defined as 
        H = (1/gamma)(a_J m x p +  b_J p)
"""
@kernel function df_torque_kernel!(@Const(m), h, gamma::T, @Const(aj), bj::T, px::T, py::T,
                                   pz::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds a = aj[id] / gamma

    b::T = bj / gamma
    @inbounds mx::T = m[j + 1]
    @inbounds my::T = m[j + 2]
    @inbounds mz::T = m[j + 3]
    @inbounds h[j + 1] = a * cross_x(mx, my, mz, px, py, pz) + b*px
    @inbounds h[j + 2] = a * cross_y(mx, my, mz, px, py, pz) + b*py
    @inbounds h[j + 3] = a * cross_z(mx, my, mz, px, py, pz) + b*pz
end

"""
The kernel torque_kernel! compute the effective field defined as 
        H = (1/gamma) m x h
"""
@kernel function torque_kernel!(@Const(m), h, gamma::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    a::T = 1/gamma

    @inbounds mx, my, mz = m[j+1], m[j+2], m[j+3]
    @inbounds hx, hy, hz = h[j+1], h[j+2], h[j+3]

    @inbounds h[j + 1] = a * cross_x(mx, my, mz, hx, hy, hz)
    @inbounds h[j + 2] = a * cross_y(mx, my, mz, hx, hy, hz)
    @inbounds h[j + 3] = a * cross_z(mx, my, mz, hx, hy, hz)
end

"""
The kernel slonczewski_torque_kernel! compute the effective field defined as 
        H = (beta*J)(epsilon* m x m_p +  xi*m_p)
"""
@kernel function slonczewski_torque_kernel!(@Const(m), h, @Const(J), lambda_sq::T, P::T,
                                            xi::T, ft::T, px::T, py::T,
                                            pz::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds mx::T = m[j + 1]
    @inbounds my::T = m[j + 2]
    @inbounds mz::T = m[j + 3]

    mp::T = mx * px + my * py + mz * pz
    epsilon::T = P * lambda_sq / (lambda_sq + 1 + (lambda_sq - 1) * mp);

    @inbounds a = J[id]*ft # note ft is multiplied by the coefficient beta 
    @inbounds h[j + 1] = a * (epsilon*cross_x(mx, my, mz, px, py, pz) + xi*px)
    @inbounds h[j + 2] = a * (epsilon*cross_y(mx, my, mz, px, py, pz) + xi*py)
    @inbounds h[j + 3] = a * (epsilon*cross_z(mx, my, mz, px, py, pz) + xi*pz)
end

"""
The kernel zhangli_torque_kernel! compute the effective field defined as 
   H = (b/gamma)*[m x (J.nabla) m + xi (J.nabla) m]
"""
@kernel function zhangli_torque_kernel!(@Const(m), h, @Const(bJ), @Const(ngbs), xi::T,
                                        ut::T, dx::T, dy::T, dz::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    fx::T, fy::T, fz::T = T(0), T(0), T(0)

    #x-direction
    i1::Int32 = ngbs[1, I] #we assume that i1<0 for the area with Ms=0
    i2::Int32 = ngbs[2, I]
    # i1 * i2 may overflow
    factor::T = (i1 > 0 && i2 > 0) ? 1 / (2 * dx) : 1 / dx
    i1 < 0 && (i1 = I)
    i2 < 0 && (i2 = I)
    j1 = 3 * i1 - 2
    j2 = 3 * i2 - 2
    @inbounds u = bJ[j] * factor
    @inbounds fx += u * (m[j2] - m[j1])
    @inbounds fy += u * (m[j2 + 1] - m[j1 + 1])
    @inbounds fz += u * (m[j2 + 2] - m[j1 + 2])

    #y-direction
    i1 = ngbs[3, I]
    i2 = ngbs[4, I]
    factor = (i1 > 0 && i2 > 0) ? 1 / (2 * dy) : 1 / dy
    i1 < 0 && (i1 = I)
    i2 < 0 && (i2 = I)
    j1 = 3 * i1 - 2
    j2 = 3 * i2 - 2
    @inbounds u = bJ[j+1] * factor
    @inbounds fx += u * (m[j2] - m[j1])
    @inbounds fy += u * (m[j2 + 1] - m[j1 + 1])
    @inbounds fz += u * (m[j2 + 2] - m[j1 + 2])

    #z-direction
    i1 = ngbs[5, I]
    i2 = ngbs[6, I]
    factor = (i1 > 0 && i2 > 0) ? 1 / (2 * dz) : 1 / dz
    i1 < 0 && (i1 = I)
    i2 < 0 && (i2 = I)
    j1 = 3 * i1 - 2
    j2 = 3 * i2 - 2
    @inbounds u = bJ[j+2] * factor
    @inbounds fx += u * (m[j2] - m[j1])
    @inbounds fy += u * (m[j2 + 1] - m[j1 + 1])
    @inbounds fz += u * (m[j2 + 2] - m[j1 + 2])

    fx = ut * fx
    fy = ut * fy
    fz = ut * fz

    # the above part is h = (b/gamma)*ut*(J.nabla) m, note we have divided ut by gamma.

    @inbounds mx::T = m[j + 0]
    @inbounds my::T = m[j + 1]
    @inbounds mz::T = m[j + 2]
    @inbounds h[j + 0] = cross_x(mx, my, mz, fx, fy, fz) + xi*fx
    @inbounds h[j + 1] = cross_y(mx, my, mz, fx, fy, fz) + xi*fy
    @inbounds h[j + 2] = cross_z(mx, my, mz, fx, fy, fz) + xi*fz
end
