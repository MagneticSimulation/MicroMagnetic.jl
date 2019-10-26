using LinearAlgebra
using CUDAnative, CuArrays

function compute_distance_kernel!(ds::CuDeviceArray{T, 1}, m1::CuDeviceArray{T, 1},
                                  m2::CuDeviceArray{T, 1},  N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        @inbounds mx1 = m1[j]
        @inbounds my1 = m1[j+1]
        @inbounds mz1 = m1[j+2]
        @inbounds mx2 = m2[j]
        @inbounds my2 = m2[j+1]
        @inbounds mz2 = m2[j+2]
        a = -my2*mz1 + my1*mz2  #m1xm2, x
        b = mx2*mz1 - mx1*mz2   #m1xm2, y
        c = -mx2*my1 + mx1*my2  #m1xm2, z
        mm = mx1*mx2 + my1*my2 + mz1*mz2
        d = CUDAnative.sqrt(a*a+b*b+c*c)
        @inbounds ds[index] = CUDAnative.atan2(d, mm)
    end
    return nothing
end

function compute_distance(neb::NEB_GPU_MPI)
    nxyz = neb.sim.nxyz
    N = neb.N
    dof = 3*nxyz
    blk, thr = CuArrays.cudims(nxyz)
    ds = neb.sim.energy #we borrow the sim.energy
    for n=0:N
       m1 = n==0 ? neb.image_l : view(neb.spin, (n-1)*dof+1:n*dof)
       m2 = n==N ? neb.image_r : view(neb.spin, n*dof+1:(n+1)*dof)
       @cuda blocks=blk threads=thr compute_distance_kernel!(ds, m1, m2, nxyz)
       neb.distance[n+1] = LinearAlgebra.norm(ds)
   end
   return nothing
end

#compute h = h - (m.h)m
function neb_projection_kernel!(h::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        @inbounds mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
        @inbounds h[j] = h[j] - mh*mx
        @inbounds h[j+1] = h[j+1] - mh*my
        @inbounds h[j+2] = h[j+2] - mh*mz
    end
    return nothing
end


function compute_tangents_kernel!(t::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                        left_m::CuDeviceArray{T, 1}, right_m::CuDeviceArray{T, 1},
                        energy::CuDeviceArray{T, 1}, N::Int64, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    dof = 3*nxyz
    if 0 < index <= dof
        for n = 1:N
            k = (n-1)*dof + index
            E1 = energy[n]
            E2 = energy[n+1]
            E3 = energy[n+2]
            dEmax = max(abs(E3-E2),abs(E2-E1))
            dEmin = min(abs(E3-E2),abs(E2-E1))

            m_kp = n==N ? right_m[index] : m[k+dof]
            m_k =  m[k]
            m_km = n==1 ? left_m[index] : m[k-dof]
            tip = m_kp - m_k
            tim = m_k - m_km

            if (E1>E2)&&(E2>E3)
              t[k] = tim
            elseif (E3>E2)&&(E2>E1)
              t[k] = tip
            elseif E3>E1
              t[k] = dEmax*tip+dEmin*tim
            elseif E3<E1
              t[k] = dEmin*tip+dEmax*tim
            else
              t[k] = tim + tip
            end
        end
    end
    return nothing
end


function compute_tangents(neb::NEB_GPU_MPI)
    N = neb.N
    nxyz = neb.sim.nxyz
    blk, thr = CuArrays.cudims(3*nxyz)
    @cuda blocks=blk threads=thr compute_tangents_kernel!(neb.tangent, neb.spin,
                                                           neb.image_l, neb.image_r,
                                                           neb.energy, N, nxyz)

    blk, thr = CuArrays.cudims(N*nxyz)
    @cuda blocks=blk threads=thr neb_projection_kernel!(neb.tangent, neb.spin, N*nxyz)

    dof = 3*nxyz
    for n=1:N
        t = view(neb.tangent, dof*(n-1)+1:dof*n)
        norm_t = LinearAlgebra.norm(t)
        t .= t/norm_t
    end

   return nothing
end

#compute LLG equation with alpha=1
function neb_llg_rhs_kernel!(dw_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                        h::CuDeviceArray{T, 1}, gamma::T, N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        @inbounds mm = mx*mx + my*my + mz*mz
        @inbounds mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
        @inbounds h1 = mm*h[j] - mh*mx
        @inbounds h2 = mm*h[j+1] - mh*my
        @inbounds h3 = mm*h[j+2] - mh*mz

        @inbounds dw_dt[j] = 0.5*gamma*h1  #we multiply a factor of 0.5 [1/(1+alpha^2) in orginal LLG equation]
        @inbounds dw_dt[j+1] = 0.5gamma*h2
        @inbounds dw_dt[j+2] = 0.5*gamma*h3
    end
    return nothing
end

function neb_llg_rhs_gpu(dw_dt::CuArray{T, 1}, m::CuArray{T, 1}, h::CuArray{T, 1},
                 gamma::T, N::Int64) where {T<:AbstractFloat}

    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr neb_llg_rhs_kernel!(dw_dt, m, h, gamma, N)
    return nothing
end
