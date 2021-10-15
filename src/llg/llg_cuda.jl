"""
CUDA kernel to compute the standard LLG equation,
    dm/dt = - gamma_L * (m x H) - alpha*gamma_L* m x (m x H)
where gamma_L = gamma/(1+alpha^2).
"""
function llg_rhs_kernel!(dm_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                        h::CuDeviceArray{T, 1}, pins::CuDeviceArray{Bool, 1}, alpha::T,
                        gamma::T, precession::Bool, N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        if pins[index]
            @inbounds dm_dt[j] = 0
            @inbounds dm_dt[j+1] = 0
            @inbounds dm_dt[j+2] = 0
            return nothing
        end

        a = -gamma/(1+alpha*alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        mm = mx*mx + my*my + mz*mz
        @inbounds mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
        @inbounds h1 = mm*h[j] - mh*mx
        @inbounds h2 = mm*h[j+1] - mh*my
        @inbounds h3 = mm*h[j+2] - mh*mz

        f1, f2, f3 = T(0), T(0), T(0)
        if precession
            f1 = cross_x(mx,my,mz, h1,h2,h3)
            f2 = cross_y(mx,my,mz, h1,h2,h3)
            f3 = cross_z(mx,my,mz, h1,h2,h3)
        end

        @inbounds dm_dt[j] = a * (f1 - h1 * alpha);
        @inbounds dm_dt[j+1] = a * (f2 - h2 * alpha);
        @inbounds dm_dt[j+2] = a * (f3 - h3 * alpha);
    end
    return nothing
end

"""
CUDA version of LLG call_back function that will be called by the integrator.
"""
function llg_call_back_gpu(sim::AbstractSimGPU, dm_dt::CuArray{T, 1}, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
    N = sim.nxyz
    driver = sim.driver
    effective_field(sim, spin, t)

    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_kernel!(dm_dt, spin, driver.field, sim.pins,
                                                 driver.alpha, driver.gamma, driver.precession, N)
    return nothing
end

"""
CUDA kernel to compute tau = (u.nabla) m.
"""
function field_stt_kernel!(m::CuDeviceArray{T, 1}, h_stt::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1},
                           ux::CuDeviceArray{T, 1}, uy::CuDeviceArray{T, 1}, uz::CuDeviceArray{T, 1}, ut::T,
                           dx::T, dy::T, dz::T, nx::Int64, ny::Int64, nz::Int64,
                           xperiodic::Bool, yperiodic::Bool, zperiodic::Bool, N::Int64)  where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if 0< index <= N
        fx, fy, fz = T(0), T(0), T(0)
        i,j,k = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])

        #x-direction
        i1 = _x_minus_one(i, index, nx, ny, nz, xperiodic, Ms)
        i2 = _x_plus_one(i, index, nx, ny, nz, xperiodic, Ms)
        factor = i1*i2>0 ? 1/(2*dx) : 1/dx
        i1 < 0 && (i1 = index)
        i2 < 0 && (i2 = index)
        j1 = 3*i1-2
        j2 = 3*i2-2
        @inbounds u = ut*ux[index]*factor
        @inbounds fx += u * (m[j2] - m[j1]);
        @inbounds fy += u * (m[j2+1] - m[j1+1]);
        @inbounds fz += u * (m[j2+2] - m[j1+2]);

        #y-direction
        i1 = _y_minus_one(j, index, nx, ny, nz, yperiodic, Ms)
        i2 = _y_plus_one(j, index, nx, ny, nz, yperiodic, Ms)
        factor = i1*i2>0 ? 1/(2*dy) : 1/dy
        i1 < 0 && (i1 = index)
        i2 < 0 && (i2 = index)
        j1 = 3*i1-2
        j2 = 3*i2-2
        @inbounds u = ut*uy[index]*factor
        @inbounds fx += u * (m[j2] - m[j1]);
        @inbounds fy += u * (m[j2+1] - m[j1+1]);
        @inbounds fz += u * (m[j2+2] - m[j1+2]);


        #z-direction
        i1 = _z_minus_one(k, index, nx, ny, nz, zperiodic, Ms)
        i2 = _z_plus_one(k, index, nx, ny, nz, zperiodic, Ms)
        factor = i1*i2>0 ? 1/(2*dz) : 1/dz
        i1 < 0 && (i1 = index)
        i2 < 0 && (i2 = index)
        j1 = 3*i1-2
        j2 = 3*i2-2
        @inbounds u = ut*uz[index]*factor
        @inbounds fx += u * (m[j2] - m[j1]);
        @inbounds fy += u * (m[j2+1] - m[j1+1]);
        @inbounds fz += u * (m[j2+2] - m[j1+2]);

        @inbounds h_stt[3*index-2] = fx
        @inbounds h_stt[3*index-1] = fy
        @inbounds h_stt[3*index] = fz

    end
    return nothing
end


"""
Wrapper for function field_stt_kernel! [to compute tau = (u.nabla) m]
"""
function compute_field_stt_gpu(m::CuArray{T, 1}, h_stt::CuArray{T, 1}, Ms::CuArray{T, 1},
                           ux::CuArray{T, 1}, uy::CuArray{T, 1}, uz::CuArray{T, 1}, ut::T,
                           dx::T, dy::T, dz::T, nx::Int64, ny::Int64, nz::Int64,
                           xperiodic::Bool, yperiodic::Bool, zperiodic::Bool, N::Int64) where {T<:AbstractFloat}

    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr field_stt_kernel!(m, h_stt, Ms,
                           ux, uy, uz, ut, dx, dy, dz, nx, ny, nz,
                           xperiodic, yperiodic, zperiodic,N)

    return nothing
end

"""
CUDA kernel to compute the STT torque and add it to dm/dt:
    dm/dt += (1+alpha*beta)/(1+alpha^2)*tau - (beta-alpha)/(1+alpha^2)*(m x tau)
where tau = (u.nabla) m, here tau is represented by h_stt.
"""
function add_stt_rhs_kernel!(dm_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                        h_stt::CuDeviceArray{T, 1}, pins::CuDeviceArray{Bool, 1}, alpha::T,
                        beta::T,  N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        if pins[index]
            return nothing
        end

        c1 = (1 + alpha* beta)/(1+alpha*alpha)
        c2 = (beta-alpha)/(1+alpha*alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        mm = mx*mx + my*my + mz*mz
        @inbounds mht = mx*h_stt[j] + my*h_stt[j+1] + mz*h_stt[j+2]
        @inbounds hp1 = mm*h_stt[j] - mht*mx;
        @inbounds hp2 = mm*h_stt[j+1] - mht*my;
        @inbounds hp3 = mm*h_stt[j+2] - mht*mz;

        mth1 = cross_x(mx,my,mz, hp1, hp2, hp3)
        mth2 = cross_y(mx,my,mz, hp1, hp2, hp3)
        mth3 = cross_z(mx,my,mz, hp1, hp2, hp3)

        @inbounds dm_dt[j] += c1 * hp1 - c2 * mth1
        @inbounds dm_dt[j+1] += c1 * hp2 - c2 * mth2
        @inbounds dm_dt[j+2] += c1 * hp3 - c2 * mth3

    end
    return nothing
end

function llg_stt_call_back_gpu(sim::AbstractSimGPU, dm_dt::CuArray{T, 1}, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
    N = sim.nxyz
    driver = sim.driver
    mesh = sim.mesh
    effective_field(sim, spin, t)

    ut = driver.ufun(t)
    compute_field_stt_gpu(spin, driver.h_stt, sim.Ms, driver.ux, driver.uy, driver.uz, T(ut),
                       mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz,
                       mesh.xperiodic, mesh.yperiodic, mesh.zperiodic, N)

    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_kernel!(dm_dt, spin, driver.field, sim.pins,
                                                 driver.alpha, driver.gamma, true, N)
    synchronize()
    @cuda blocks=blk threads=thr add_stt_rhs_kernel!(dm_dt, spin, driver.h_stt, sim.pins, driver.alpha, driver.beta, N)

    return nothing

end

"""
CUDA kernel to compute the STT torque for CPP case and add it to dm/dt:
    dm/dt += (a_J+alpha*b_J)/(1+alpha^2)*p_perp - (b_J-alpha*a_J)/(1+alpha^2)*(m x p_perp)
where p_perp = p - (m.p)m
"""
function add_cpp_rhs_kernel!(dm_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, pins::CuDeviceArray{Bool, 1},
                        alpha::T, a_J::CuDeviceArray{T, 1}, bj::T, px::T, py::T, pz::T,
                        N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        if pins[index]
            return nothing
        end

        @inbounds aj = a_J[index]
        c1 = (aj + alpha*bj)/(1+alpha*alpha)
        c2 = (bj - alpha*aj)/(1+alpha*alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        mm = mx*mx + my*my + mz*mz
        mpt = mx*px + my*py + mz*pz
        @inbounds hp1 = mm*px-mpt*mx;
        @inbounds hp2 = mm*py-mpt*my;
        @inbounds hp3 = mm*pz-mpt*mz;

        mth1 = cross_x(mx, my, mz, hp1, hp2, hp3)
        mth2 = cross_y(mx, my, mz, hp1, hp2, hp3)
        mth3 = cross_z(mx, my, mz, hp1, hp2, hp3)

        @inbounds dm_dt[j] += c1 * hp1 - c2 * mth1
        @inbounds dm_dt[j+1] += c1 * hp2 - c2 * mth2
        @inbounds dm_dt[j+2] += c1 * hp3 - c2 * mth3

    end
    return nothing
end


function llg_cpp_call_back_gpu(sim::AbstractSimGPU, dm_dt::CuArray{T, 1}, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
    N = sim.nxyz
    driver = sim.driver
    mesh = sim.mesh
    effective_field(sim, spin, t)

    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_kernel!(dm_dt, spin, driver.field, sim.pins,
                                                 driver.alpha, driver.gamma, true, N)
    synchronize()
    px, py, pz = T(driver.p[1]), T(driver.p[2]), T(driver.p[3])
    @cuda blocks=blk threads=thr add_cpp_rhs_kernel!(dm_dt, spin, sim.pins, driver.alpha, driver.aj, driver.bj, px, py, pz, N)

    return nothing

end

function llg_stt_cpp_call_back_gpu(sim::AbstractSimGPU, dm_dt::CuArray{T, 1}, spin::CuArray{T, 1}, t::Float64) where {T<:AbstractFloat}
    N = sim.nxyz
    driver = sim.driver
    mesh = sim.mesh
    effective_field(sim, spin, t)

    ut = driver.ufun(t)
    compute_field_stt_gpu(spin, driver.h_stt, sim.Ms, driver.ux, driver.uy, driver.uz, T(ut),
                       mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz,
                       mesh.xperiodic, mesh.yperiodic, mesh.zperiodic, N)

    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_kernel!(dm_dt, spin, driver.field, sim.pins,
                                                 driver.alpha, driver.gamma, true, N)
    synchronize()
    @cuda blocks=blk threads=thr add_stt_rhs_kernel!(dm_dt, spin, driver.h_stt, sim.pins, driver.alpha, driver.beta, N)

    synchronize()
    px, py, pz = T(driver.p[1]), T(driver.p[2]), T(driver.p[3])
    @cuda blocks=blk threads=thr add_cpp_rhs_kernel!(dm_dt, spin, sim.pins, driver.alpha, driver.aj, driver.bj, px, py, pz, N)

    return nothing

end
