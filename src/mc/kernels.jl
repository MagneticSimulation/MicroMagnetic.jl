@kernel function uniform_random_sphere_kernel!(m, @Const(rnd))
    I = @index(Global)
    j = 3 * I - 2
    @inbounds phi = rnd[j] * 2 * pi
    @inbounds ct = 2 * rnd[j + 1] - 1
    st = sqrt(1 - ct * ct)
    @inbounds m[j] = st * cos(phi)
    @inbounds m[j + 1] = st * sin(phi)
    @inbounds m[j + 2] = ct
end

@kernel function uniform_random_circle_xy_kernel!(m, @Const(rnd))
    I = @index(Global)
    j = 3 * I - 2
    @inbounds phi = rnd[j] * 2 * pi
    @inbounds m[j] = cos(phi)
    @inbounds m[j + 1] = sin(phi)
    @inbounds m[j + 2] = 0
end

@kernel function dE_zeeman_anisotropy_energy_kernel!(@Const(m), @Const(next_m),
                                                     @Const(shape), energy, Hx::T, Hy::T,
                                                     Hz::T, Ku::T, Kc::T, ux::T, uy::T,
                                                     uz::T, bias::Int64,
                                                     cubic::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    a, b, c = @index(Global, NTuple)
    @inbounds in_shape = shape[I]
    if in_shape
        sign = cubic ? 1 : -1
        if mod(a + sign * b + c, 3) == bias
            i = 3 * I - 2
            @inbounds dmx = next_m[i] - m[i]
            @inbounds dmy = next_m[i + 1] - m[i + 1]
            @inbounds dmz = next_m[i + 2] - m[i + 2]

            delta_E = -(dmx * Hx + dmy * Hy + dmz * Hz) #zeeman

            @inbounds delta_E += Ku * (m[i] * ux + m[i + 1] * uy + m[i + 2] * uz)^2 #Anisotropy
            @inbounds delta_E -= Ku *
                                 (next_m[i] * ux + next_m[i + 1] * uy + next_m[i + 2] * uz)^2
            @inbounds delta_E += Kc * (m[i + 2]^4 - next_m[i + 2]^4) #Cubic Anisotropy for z
            @inbounds delta_E += Kc * (m[i + 1]^4 - next_m[i + 1]^4) #Cubic Anisotropy for y
            @inbounds delta_E += Kc * (m[i]^4 - next_m[i]^4) #Cubic Anisotropy for x
            @inbounds energy[I] = delta_E
        end
    end
end

@kernel function zeeman_anisotropy_energy_kernel!(@Const(m), @Const(shape), energy, Hx::T,
                                                  Hy::T, Hz::T, Ku::T, Kc::T, ux::T, uy::T,
                                                  uz::T) where {T<:AbstractFloat}
    I = @index(Global)
    @inbounds in_shape = shape[I]
    if in_shape
        i = 3 * I - 2
        @inbounds mx = m[i]
        @inbounds my = m[i + 1]
        @inbounds mz = m[i + 2]

        delta_E::T = 0
        delta_E = -(mx * Hx + my * Hy + mz * Hz) #zeeman
        delta_E -= Ku * (mx * ux + my * uy + mz * uz)^2 #Anisotropy
        delta_E += Kc * (mx^4 + my^4 + mz^4)
        @inbounds energy[I] = delta_E
    end
end

@kernel function add_dE_exch_dmi_energy_kernel!(@Const(m), @Const(next_m), @Const(shape),
                                                energy, @Const(ngbs), n_ngbs::Int64, Jx::T,
                                                Jy::T, Jz::T, @Const(D0), bias::Int64,
                                                cubic::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    a, b, c = @index(Global, NTuple)
    @inbounds in_shape = shape[I]
    if in_shape
        sign = cubic ? 1 : -1
        if mod(a + sign * b + c, 3) == bias
            i = 3 * I - 2
            @inbounds dmx = next_m[i] - m[i]
            @inbounds dmy = next_m[i + 1] - m[i + 1]
            @inbounds dmz = next_m[i + 2] - m[i + 2]

            delta_E::T = 0
            for j in 1:n_ngbs
                @inbounds id = ngbs[j, I]

                @inbounds Dx = D0[1, j]
                @inbounds Dy = D0[2, j]
                @inbounds Dz = D0[3, j]

                if id > 0 && shape[id]
                    k = 3 * id - 2
                    @inbounds sx = m[k]
                    @inbounds sy = m[k + 1]
                    @inbounds sz = m[k + 2]
                    delta_E -= (Jx * dmx * sx + Jy * dmy * sy + Jz * dmz * sz) #exchange
                    delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
                end
            end
            @inbounds energy[I] += delta_E
        end
    end
end

@kernel function add_exch_dmi_energy_kernel!(@Const(m), @Const(shape),
                                             energy, @Const(ngbs), n_ngbs::Int64, Jx::T,
                                             Jy::T, Jz::T,
                                             @Const(D0)) where {T<:AbstractFloat}
    I = @index(Global)
    @inbounds in_shape = shape[I]
    if in_shape
        i = 3 * I - 2
        @inbounds mx = m[i]
        @inbounds my = m[i + 1]
        @inbounds mz = m[i + 2]

        delta_E::T = 0
        for j in 1:n_ngbs
            @inbounds id = ngbs[j, I]

            @inbounds Dx = D0[1, j]
            @inbounds Dy = D0[2, j]
            @inbounds Dz = D0[3, j]

            if id > 0 && shape[id]
                k = 3 * id - 2
                @inbounds sx = m[k]
                @inbounds sy = m[k + 1]
                @inbounds sz = m[k + 2]
                delta_E -= (Jx * mx * sx + Jy * my * sy + Jz * mz * sz) #exchange
                delta_E += volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end
        end
        @inbounds energy[I] += delta_E
    end
end

@kernel function run_monte_carlo_kernel!(m, @Const(next_m), @Const(rnd), @Const(shape),
                                         @Const(delta_enengy), temp::T, bias::Int64,
                                         cubic::Bool) where {T<:AbstractFloat}
    I = @index(Global)
    a, b, c = @index(Global, NTuple)
    @inbounds in_shape = shape[I]
    if in_shape
        sign = cubic ? 1 : -1
        # bias should be 0, 1 and 2
        # We only deal with the specfic sites 
        if mod(a + sign * b + c, 3) == bias
            i = 3 * I - 2
            @inbounds rnd_i = rnd[i + 2]

            delta_E::T = 0
            @inbounds delta_E = delta_enengy[I]
            if delta_E < 0 || rnd_i < exp(-delta_E / temp)
                @inbounds m[i] = next_m[i]
                @inbounds m[i + 1] = next_m[i + 1]
                @inbounds m[i + 2] = next_m[i + 2]
            end
        end
    end
end
