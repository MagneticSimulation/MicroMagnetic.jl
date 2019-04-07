function llg_call_back_gpu(sim::MicroSimGPU, dw_dt::CuArray{T, 1}, t::Float64, omega::CuArray{T, 1}) where {T<:AbstractFloat}

  #dw_dt = sim.driver.ode.dw_dt
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)
  llg_rhs_gpu(dw_dt, sim.spin, sim.driver.field, omega, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.nxyz)
  return nothing

end


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
