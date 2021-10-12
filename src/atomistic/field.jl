function effective_field(zeeman::ZeemanGPU, sim::AtomicSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

    function __kernel!(m, h, h_static, energy, mu_s, N)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        if 0 < i <= N
            j = 3*(i-1)
            @inbounds h[j+1] = h_static[j+1]
            @inbounds h[j+2] = h_static[j+2]
            @inbounds h[j+3] = h_static[j+3]
            @inbounds energy[i] = -mu_s[i]*(m[j+1]*h[j+1] + m[j+2]*h[j+2] + m[j+3]*h[j+3])
        end
        return nothing
    end

    N = sim.nxyz
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr __kernel!(spin, sim.field, zeeman.cufield, sim.energy, sim.mu_s, N)

    return nothing
end

function effective_field(anis::AnisotropyGPU, sim::AtomicSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

    function __kernel!(m, h, energy, ku, ax, ay, az, mu_s, N)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        if 0 < i <= N
            j = 3*(i-1)

            @inbounds ms_local = mu_s[i]
            if ms_local == 0.0
                @inbounds energy[i] = 0
                @inbounds h[j+1] = 0
                @inbounds h[j+2] = 0
                @inbounds h[j+3] = 0
                return nothing
            end
            ms_inv = 2.0/ms_local
            @inbounds sa = m[j+1]*ax+m[j+2]*ay+m[j+3]*az
            @inbounds h[j+1] = ku[i]*sa*ax*ms_inv
            @inbounds h[j+2] = ku[i]*sa*ay*ms_inv
            @inbounds h[j+3] = ku[i]*sa*az*ms_inv
            @inbounds energy[i] = -ku[i]*sa*sa
        end

        return nothing
    end

    N = sim.nxyz
    axis = anis.axis
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr __kernel!(spin, sim.field, sim.energy, anis.Ku, axis[1], axis[2], axis[3], sim.mu_s, N)

    return nothing
end


function effective_field(exch::HeisenbergExchange, sim::AtomicSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

    function __kernel!(m, h, energy, Js, ngbs, n_ngbs, mu_s, N)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        if 0 < i <= N
            j = 3*(i-1)
            @inbounds ms_local = mu_s[i]
            if ms_local == 0.0
                @inbounds energy[i] = 0
                @inbounds h[j+1] = 0
                @inbounds h[j+2] = 0
                @inbounds h[j+3] = 0
                return nothing
            end
            ms_inv = T(1/ms_local)

            fx, fy, fz = T(0), T(0), T(0)
            for k = 1:n_ngbs
                @inbounds id = ngbs[k, i]
                if id>0 && mu_s[id] > 0
                    x = 3*id-2
                    @inbounds fx += Js[k]*m[x]
                    @inbounds fy += Js[k]*m[x+1]
                    @inbounds fz += Js[k]*m[x+2]
                end
            end
            @inbounds energy[i] = -0.5*(fx*m[j+1] + fy*m[j+2] + fz*m[j+3])
            @inbounds h[j+1] = fx*ms_inv
            @inbounds h[j+2] = fy*ms_inv
            @inbounds h[j+3] = fz*ms_inv
        end

        return nothing
    end

    N = sim.nxyz
    blk, thr = cudims(N)
    ngbs = sim.mesh.ngbs
    n_ngbs = sim.mesh.n_ngbs
    @cuda blocks=blk threads=thr __kernel!(spin, sim.field, sim.energy, exch.Js, ngbs, n_ngbs, sim.mu_s, N)

    return nothing
end


function effective_field(dmi::HeisenbergBulkDMI, sim::AtomicSimGPU, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}

    function __kernel!(m, h, energy, D, ngbs, n_ngbs, mu_s, N)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        
        ax = (1.0,-1.0, 0.0, 0.0, 0.0, 0.0)
        ay = (0.0, 0.0, 1.0,-1.0, 0.0, 0.0)
        az = (0.0, 0.0, 0.0, 0.0, 1.0,-1.0)

        if 0 < i <= N
            j = 3*(i-1)
            @inbounds ms_local = mu_s[i]
            if ms_local == 0.0
                @inbounds energy[i] = 0
                @inbounds h[j+1] = 0
                @inbounds h[j+2] = 0
                @inbounds h[j+3] = 0
                return nothing
            end
            ms_inv = T(1/ms_local)

            fx, fy, fz = T(0), T(0), T(0)

            for k = 1:n_ngbs
                @inbounds id = ngbs[k, i]
                if id>0 && mu_s[id] > 0
                    x = 3*id-2
                    @inbounds fx += D*cross_x(ax[k],ay[k],az[k],m[x],m[x+1],m[x+2]);
                    @inbounds fy += D*cross_y(ax[k],ay[k],az[k],m[x],m[x+1],m[x+2]);
                    @inbounds fz += D*cross_z(ax[k],ay[k],az[k],m[x],m[x+1],m[x+2]);
                end
            end
            @inbounds energy[i] = -0.5*(fx*m[j+1] + fy*m[j+2] + fz*m[j+3])
            @inbounds h[j+1] = fx*ms_inv
            @inbounds h[j+2] = fy*ms_inv
            @inbounds h[j+3] = fz*ms_inv
        end

        return nothing
    end

    N = sim.nxyz
    blk, thr = cudims(N)
    ngbs = sim.mesh.ngbs
    n_ngbs = sim.mesh.n_ngbs
    @cuda blocks=blk threads=thr __kernel!(spin, sim.field, sim.energy, dmi.D, ngbs, n_ngbs, sim.mu_s, N)

    return nothing
end
