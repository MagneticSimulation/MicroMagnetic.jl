
@doc raw"""
compute the exchange interaction h_ex,  which is defined as 

```math
\mathbf{H}_\mathrm{ex, i} = \frac{1}{\mu_s} \sum_{\langle i, j\rangle} J_{i,j} \mathbf{m}_{j}.
```

Jij is represented by Js, which is an array with its length equals to n_ngbs

Note that the return value is h * alpha + h_ex 
"""
@kernel function atomistic_exchange_kernel!(h, energy, @Const(Js), @Const(m), @Const(mu_s),
                                            @Const(ngbs), n_ngbs,
                                            alpha::T) where {T<:AbstractFloat}
    I = @index(Global)
    i = 3 * I - 2

    @inbounds ms_local = mu_s[I]
    if ms_local == 0.0
        @inbounds energy[I] = 0
        @inbounds h[i + 0] = 0
        @inbounds h[i + 1] = 0
        @inbounds h[i + 2] = 0
    else
        ms_inv::T = 1 / ms_local

        fx, fy, fz = T(0), T(0), T(0)
        for j in 1:n_ngbs
            @inbounds id = ngbs[j, I]
            if id > 0 && mu_s[id] > 0
                k = 3 * id - 2
                @inbounds fx += (Js[j]) * m[k + 0]
                @inbounds fy += (Js[j]) * m[k + 1]
                @inbounds fz += (Js[j]) * m[k + 2]
            end
        end

        @inbounds energy_I::T = -0.5 * (fx * m[i] + fy * m[i + 1] + fz * m[i + 2])
        @inbounds energy[I] = energy[I] * alpha + energy_I
        @inbounds h[i + 0] = h[i + 0] * alpha + fx * ms_inv
        @inbounds h[i + 1] = h[i + 1] * alpha + fy * ms_inv
        @inbounds h[i + 2] = h[i + 2] * alpha + fz * ms_inv
    end
end

@doc raw"""
compute the dmi interaction h_dmi,  which is defined as 

```math
\mathcal{H}_\mathrm{dmi} =  \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right)
```

Jij is a 2d array with size (3, n_ngbs)
"""
@kernel function atomistic_dmi_kernel!(h::AbstractArray{T,1}, energy::AbstractArray{T,1},
                                       @Const(Dij), @Const(m), @Const(mu_s), @Const(ngbs),
                                       n_ngbs) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds ms_local::T = mu_s[I]
    if ms_local == 0.0
        @inbounds energy[I] = 0
        @inbounds h[j] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
    else
        ms_inv::T = 1 / ms_local

        fx, fy, fz = T(0), T(0), T(0)
        for k in 1:n_ngbs
            @inbounds id = ngbs[k, I]
            @inbounds if id > 0 && mu_s[id] > 0
                x = 3 * id - 2
                @inbounds Dx::T = Dij[1, k]
                @inbounds Dy::T = Dij[2, k]
                @inbounds Dz::T = Dij[3, k]
                @inbounds fx += cross_x(Dx, Dy, Dz, m[x], m[x + 1], m[x + 2])
                @inbounds fy += cross_y(Dx, Dy, Dz, m[x], m[x + 1], m[x + 2])
                @inbounds fz += cross_z(Dx, Dy, Dz, m[x], m[x + 1], m[x + 2])
            end
        end
        @inbounds energy[I] = -0.5 * (fx * m[j] + fy * m[j + 1] + fz * m[j + 2])
        @inbounds h[j] = fx * ms_inv
        @inbounds h[j + 1] = fy * ms_inv
        @inbounds h[j + 2] = fz * ms_inv
    end
end


@doc raw"""
In canted AFM DMI has an pmpm pattern
compute the dmi interaction h_dmi,  which is defined as 

```math
\mathcal{H}_{\mathrm{dmi}}=\sum_{i}\left(-1\right)^{i}\mathbf{D}_{i}\cdot\left(\mathbf{m}_{i}\times\mathbf{m}_{i+1}\right)
```

Jij is a 2d array with size (3, n_ngbs)
"""
@kernel function atomistic_canted_dmi_kernel!(h::AbstractArray{T,1}, energy::AbstractArray{T,1},
    @Const(Dij), @Const(m), @Const(mu_s), @Const(ngbs),
    n_ngbs) where {T<:AbstractFloat}
    I = @index(Global)
    ii, jj, kk = @index(Global, NTuple)
    j = 3 * I - 2

    pmf=mod(ii+jj+kk,2)*2-1

    @inbounds ms_local::T = mu_s[I]
    if ms_local == 0.0
        @inbounds energy[I] = 0
        @inbounds h[j] = 0
        @inbounds h[j+1] = 0
        @inbounds h[j+2] = 0
    else
        ms_inv::T = 1 / ms_local

        fx, fy, fz = T(0), T(0), T(0)
        for k in 1:n_ngbs
            @inbounds id = ngbs[k, I]
            @inbounds if id > 0 && mu_s[id] > 0
                x = 3 * id - 2
                @inbounds Dx::T = pmf*Dij[1, k]
                @inbounds Dy::T = pmf*Dij[2, k]
                @inbounds Dz::T = pmf*Dij[3, k]
                @inbounds fx += cross_x(Dx, Dy, Dz, m[x], m[x+1], m[x+2])
                @inbounds fy += cross_y(Dx, Dy, Dz, m[x], m[x+1], m[x+2])
                @inbounds fz += cross_z(Dx, Dy, Dz, m[x], m[x+1], m[x+2])
            end
        end
        @inbounds energy[I] = -0.5 * (fx * m[j] + fy * m[j+1] + fz * m[j+2])
        @inbounds h[j] = fx * ms_inv
        @inbounds h[j+1] = fy * ms_inv
        @inbounds h[j+2] = fz * ms_inv
    end
end


@kernel function tube_bulk_dmi_kernel!(h, energy, D::T, @Const(Dij), @Const(m),
                                       @Const(mu_s), @Const(ngbs), n_ngbs,
                                       nr) where {T<:AbstractFloat}
    i = @index(Global)
    j = 3 * (I - 1)
    @inbounds ms_local::T = mu_s[i]
    if ms_local == 0.0
        @inbounds energy[I] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        ms_inv::T = 1 / ms_local

        fx, fy, fz = T(0), T(0), T(0)

        # the order of Dij is r_{n_r, 1}, r12, r23, ..., r_{(n_r-1)n_r},  r_{n_r, 1}
        # the left neighbour
        @inbounds id = ngbs[1, i]
        if id > 0 && mu_s[id] > 0
            x = 3 * id - 2
            k = i % nr  #should be i rather than id
            k == 0 && (k += nr) # the rangle of k is [1, nr]
            @inbounds fx -= cross_x(Dij[1, k], Dij[2, k], Dij[3, k], m[x], m[x + 1],
                                    m[x + 2])
            @inbounds fy -= cross_y(Dij[1, k], Dij[2, k], Dij[3, k], m[x], m[x + 1],
                                    m[x + 2])
            @inbounds fz -= cross_z(Dij[1, k], Dij[2, k], Dij[3, k], m[x], m[x + 1],
                                    m[x + 2])
        end

        # the right neighbour
        @inbounds id = ngbs[2, i]
        if id > 0 && mu_s[id] > 0
            x = 3 * id - 2
            k = i % nr
            k == 0 && (k += nr)
            k += 1  # the range of k is [2, nr+1]
            @inbounds fx += cross_x(Dij[1, k], Dij[2, k], Dij[3, k], m[x], m[x + 1],
                                    m[x + 2])
            @inbounds fy += cross_y(Dij[1, k], Dij[2, k], Dij[3, k], m[x], m[x + 1],
                                    m[x + 2])
            @inbounds fz += cross_z(Dij[1, k], Dij[2, k], Dij[3, k], m[x], m[x + 1],
                                    m[x + 2])
        end

        # the bottom neighbour
        @inbounds id = ngbs[3, i]
        if id > 0 && mu_s[id] > 0
            x = 3 * id - 2
            @inbounds fx += cross_x(T(0), T(0), -D, m[x], m[x + 1], m[x + 2])
            @inbounds fy += cross_y(T(0), T(0), -D, m[x], m[x + 1], m[x + 2])
            @inbounds fz += cross_z(T(0), T(0), -D, m[x], m[x + 1], m[x + 2])
        end

        # the top neighbour
        @inbounds id = ngbs[4, i]
        if id > 0 && mu_s[id] > 0
            x = 3 * id - 2
            @inbounds fx += cross_x(T(0), T(0), D, m[x], m[x + 1], m[x + 2])
            @inbounds fy += cross_y(T(0), T(0), D, m[x], m[x + 1], m[x + 2])
            @inbounds fz += cross_z(T(0), T(0), D, m[x], m[x + 1], m[x + 2])
        end

        @inbounds energy[i] = -0.5 * (fx * m[j + 1] + fy * m[j + 2] + fz * m[j + 3])
        @inbounds h[j + 1] = fx * ms_inv
        @inbounds h[j + 2] = fy * ms_inv
        @inbounds h[j + 3] = fz * ms_inv
    end
end
