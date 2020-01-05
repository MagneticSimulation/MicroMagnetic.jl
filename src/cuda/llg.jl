function llg_rhs_kernal!(dw_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                        h::CuDeviceArray{T, 1}, omega::CuDeviceArray{T, 1}, alpha::T,
                        gamma::T, precession::Bool, N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        a = -gamma/(1+alpha*alpha)
        b = alpha*a
        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        @inbounds mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
        @inbounds h1 = h[j] - mh*mx
        @inbounds h2 = h[j+1] - mh*my
        @inbounds h3 = h[j+2] - mh*mz
        f1 = -a*h1*precession - b*cross_x(mx,my,mz, h1,h2,h3)
        f2 = -a*h2*precession - b*cross_y(mx,my,mz, h1,h2,h3)
        f3 = -a*h3*precession - b*cross_z(mx,my,mz, h1,h2,h3)

        @inbounds wx = omega[j]
        @inbounds wy = omega[j+1]
        @inbounds wz = omega[j+2]

        wf = wx*f1 + wy*f2 + wz*f3
        @inbounds dw_dt[j] = f1 - 0.5*cross_x(wx, wy, wz, f1, f2, f3) + 0.25*wf*wx
        @inbounds dw_dt[j+1] = f2 - 0.5*cross_y(wx, wy, wz, f1, f2, f3) + 0.25*wf*wy
        @inbounds dw_dt[j+2] = f3 - 0.5*cross_z(wx, wy, wz, f1, f2, f3) + 0.25*wf*wz
    end
    return nothing
end

function llg_rhs_gpu(dw_dt::CuArray{T, 1}, m::CuArray{T, 1}, h::CuArray{T, 1},
                 omega::CuArray{T, 1}, alpha::T, gamma::T, precession::Bool, N::Int64) where {T<:AbstractFloat}

    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_kernal!(dw_dt, m, h, omega, alpha, gamma, precession, N)
    return nothing
end

function field_stt_kernal!(m::CuDeviceArray{T, 1}, h_stt::CuDeviceArray{T, 1}, Ms::CuDeviceArray{T, 1},
                           ux::CuDeviceArray{T, 1}, uy::CuDeviceArray{T, 1}, uz::CuDeviceArray{T, 1}, ut::T,
                           dx::T, dy::T, dz::T, nx::Int64, ny::Int64, nz::Int64,
                           xperiodic::Bool, yperiodic::Bool, zperiodic::Bool, N::Int64)  where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if 0< index <= N
        fx, fy, fz = T(0), T(0), T(0)
        i,j,k = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])

        #x-direction
        i1 = _x_minus_one(i, index, nx, ny, nz, xperiodic, Ms)
        i2 = _x_plus_one(i, index, nx, ny, nz, xperiodic, Ms)
        factor = i1*i2>0 ? 1/(2*dx) : 1/dx
        i1 < 0 && (i1 = i)
        i2 < 0 && (i2 = i)
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
        i1 < 0 && (i1 = i)
        i2 < 0 && (i2 = i)
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
        i1 < 0 && (i1 = i)
        i2 < 0 && (i2 = i)
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


#compute (\vec{u} \cdot \nabla \vec{m})
function compute_field_stt_gpu(m::CuArray{T, 1}, h_stt::CuArray{T, 1}, Ms::CuArray{T, 1},
                           ux::CuArray{T, 1}, uy::CuArray{T, 1}, uz::CuArray{T, 1}, ut::T,
                           dx::T, dy::T, dz::T, nx::Int64, ny::Int64, nz::Int64,
                           xperiodic::Bool, yperiodic::Bool, zperiodic::Bool, N::Int64) where {T<:AbstractFloat}

    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr field_stt_kernal!(m, h_stt, Ms,
                           ux, uy, uz, ut, dx, dy, dz, nx, ny, nz,
                           xperiodic, yperiodic, zperiodic,N)

    return nothing
end

function llg_rhs_stt_kernal!(dw_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                     h::CuDeviceArray{T, 1}, h_stt::CuDeviceArray{T, 1},
                     omega::CuDeviceArray{T, 1}, alpha::T, beta::T,
                     gamma::T, N::Int64) where { T<:AbstractFloat }
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        a = gamma/(1+alpha*alpha)
        b = alpha*a
        u::T = 1.0/(1+alpha*alpha)

        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        @inbounds mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
        @inbounds h1 = h[j] - mh*mx
        @inbounds h2 = h[j+1] - mh*my
        @inbounds h3 = h[j+2] - mh*mz

        @inbounds mht = mx*h_stt[j] + my*h_stt[j+1] + mz*h_stt[j+2]
        @inbounds ht1 = h_stt[j] - mht*mx;
        @inbounds ht2 = h_stt[j+1] - mht*my;
        @inbounds ht3 = h_stt[j+2] - mht*mz;

        f1 = a*h1 + u*(beta-alpha)*ht1
	    f2 = a*h2 + u*(beta-alpha)*ht2
	    f3 = a*h3 + u*(beta-alpha)*ht3

        hpt1 = b*h1 + u*(1+alpha*beta)*ht1
        hpt2 = b*h2 + u*(1+alpha*beta)*ht2
        hpt3 = b*h3 + u*(1+alpha*beta)*ht3

        f1 += cross_x(mx,my,mz, hpt1, hpt2, hpt3)
        f2 += cross_y(mx,my,mz, hpt1, hpt2, hpt3)
        f3 += cross_z(mx,my,mz, hpt1, hpt2, hpt3)

        @inbounds wx = omega[j]
        @inbounds wy = omega[j+1]
        @inbounds wz = omega[j+2]

        wf = wx*f1 + wy*f2 + wz*f3
        @inbounds dw_dt[j] = f1 - 0.5*cross_x(wx, wy, wz, f1, f2, f3) + 0.25*wf*wx
        @inbounds dw_dt[j+1] = f2 - 0.5*cross_y(wx, wy, wz, f1, f2, f3) + 0.25*wf*wy
        @inbounds dw_dt[j+2] = f3 - 0.5*cross_z(wx, wy, wz, f1, f2, f3) + 0.25*wf*wz
    end
    return nothing
end

function llg_rhs_stt_gpu(dw_dt::CuArray{T, 1}, m::CuArray{T, 1}, h::CuArray{T, 1}, h_stt::CuArray{T, 1},
                 omega::CuArray{T, 1}, alpha::T, beta::T, gamma::T, N::Int64) where {T<:AbstractFloat}

    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_stt_kernal!(dw_dt, m, h, h_stt, omega, alpha, beta, gamma, N)
    return nothing
end

function llg_rhs_stt_cpp_kernal!(dw_dt::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1},
                     h::CuDeviceArray{T, 1}, h_stt::CuDeviceArray{T, 1},
                     omega::CuDeviceArray{T, 1}, aj::CuDeviceArray{T, 1},
                     px::T, py::T, pz::T, alpha::T, beta::T, bj::T,
                     gamma::T, N::Int64) where { T<:AbstractFloat }
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index-2
        a = gamma/(1+alpha*alpha)
        b = alpha*a
        u::T = 1.0/(1+alpha*alpha)

        @inbounds aj_ = aj[index]

        @inbounds mx = m[j]
        @inbounds my = m[j+1]
        @inbounds mz = m[j+2]
        @inbounds mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
        @inbounds h1 = h[j] - mh*mx
        @inbounds h2 = h[j+1] - mh*my
        @inbounds h3 = h[j+2] - mh*mz

        @inbounds mht = mx*h_stt[j] + my*h_stt[j+1] + mz*h_stt[j+2]
        @inbounds ht1 = h_stt[j] - mht*mx;
        @inbounds ht2 = h_stt[j+1] - mht*my;
        @inbounds ht3 = h_stt[j+2] - mht*mz;

        mpt = mx*px + my*py + mz*pz
        @inbounds hp1 = px-mpt*mx;
        @inbounds hp2 = py-mpt*my;
        @inbounds hp3 = pz-mpt*mz;

        f1 = a*h1 + u*(beta-alpha)*ht1 + u*(bj-alpha*aj_)*hp1
        f2 = a*h2 + u*(beta-alpha)*ht2 + u*(bj-alpha*aj_)*hp2
        f3 = a*h3 + u*(beta-alpha)*ht3 + u*(bj-alpha*aj_)*hp3

        hpt1 = b*h1 + u*(1+alpha*beta)*ht1 + u*(aj_+alpha*bj)*hp1
        hpt2 = b*h2 + u*(1+alpha*beta)*ht2 + u*(aj_+alpha*bj)*hp2
        hpt3 = b*h3 + u*(1+alpha*beta)*ht3 + u*(aj_+alpha*bj)*hp3

        f1 += cross_x(mx,my,mz, hpt1, hpt2, hpt3)
        f2 += cross_y(mx,my,mz, hpt1, hpt2, hpt3)
        f3 += cross_z(mx,my,mz, hpt1, hpt2, hpt3)

        @inbounds wx = omega[j]
        @inbounds wy = omega[j+1]
        @inbounds wz = omega[j+2]

        wf = wx*f1 + wy*f2 + wz*f3
        @inbounds dw_dt[j] = f1 - 0.5*cross_x(wx, wy, wz, f1, f2, f3) + 0.25*wf*wx
        @inbounds dw_dt[j+1] = f2 - 0.5*cross_y(wx, wy, wz, f1, f2, f3) + 0.25*wf*wy
        @inbounds dw_dt[j+2] = f3 - 0.5*cross_z(wx, wy, wz, f1, f2, f3) + 0.25*wf*wz
    end
    return nothing
end

function llg_rhs_stt_cpp_gpu(dw_dt::CuArray{T, 1}, m::CuArray{T, 1}, h::CuArray{T, 1}, h_stt::CuArray{T, 1},
                 omega::CuArray{T, 1}, aj::CuArray{T, 1}, p::Tuple, alpha::T, beta::T, bj::T, gamma::T, N::Int64) where {T<:AbstractFloat}
    blk, thr = CuArrays.cudims(N)
    px, py, pz = T(p[1]), T(p[2]), T(p[3])
    @cuda blocks=blk threads=thr llg_rhs_stt_cpp_kernal!(dw_dt, m, h, h_stt, omega, aj, px, py, pz, alpha, beta, bj, gamma, N)
    return nothing
end

function llg_call_back_gpu(sim::MicroSimGPU, dw_dt::CuArray{T, 1}, t::Float64, omega::CuArray{T, 1}) where {T<:AbstractFloat}
    omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
    effective_field(sim, sim.spin, t)
    llg_rhs_gpu(dw_dt, sim.spin, sim.driver.field, omega, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.nxyz)
    return nothing
end

function llg_stt_call_back_gpu(sim::AbstractSim, dw_dt::CuArray{T, 1}, t::Float64, omega::CuArray{T, 1}) where {T<:AbstractFloat}

    driver = sim.driver
    mesh = sim.mesh
    omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
    effective_field(sim, sim.spin, t)
    ut = driver.ufun(t)
    compute_field_stt_gpu(sim.spin, driver.h_stt, sim.Ms, driver.ux, driver.uy, driver.uz, T(ut),
                       mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz,
                       mesh.xperiodic, mesh.yperiodic, mesh.zperiodic, sim.nxyz)
    llg_rhs_stt_gpu(dw_dt, sim.spin, driver.field, driver.h_stt, omega, driver.alpha, driver.beta, driver.gamma, sim.nxyz)

    return nothing

end

function llg_stt_cpp_call_back_gpu(sim::AbstractSim, dw_dt::CuArray{T, 1}, t::Float64, omega::CuArray{T, 1}) where {T<:AbstractFloat}

    driver = sim.driver
    mesh = sim.mesh
    omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
    effective_field(sim, sim.spin, t)
    ut = driver.ufun(t)
    compute_field_stt_gpu(sim.spin, driver.h_stt, sim.Ms, driver.ux, driver.uy, driver.uz, T(ut),
                       mesh.dx, mesh.dy, mesh.dz, mesh.nx, mesh.ny, mesh.nz,
                       mesh.xperiodic, mesh.yperiodic, mesh.zperiodic, sim.nxyz)

    llg_rhs_stt_cpp_gpu(dw_dt, sim.spin, driver.field, driver.h_stt, omega, driver.aj, driver.p, driver.alpha,
                        driver.beta, driver.bj, driver.gamma, sim.nxyz)

    return nothing

end
