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
        mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
        h1 = h[j] - mh*m[j]
        h2 = h[j+1] - mh*m[j+1]
        h3 = h[j+2] - mh*m[j+2]
        f1 = -a*h1*precession - b*cross_x(m[j],m[j+1],m[j+2], h1,h2,h3)
        f2 = -a*h2*precession - b*cross_y(m[j],m[j+1],m[j+2], h1,h2,h3)
        f3 = -a*h3*precession - b*cross_z(m[j],m[j+1],m[j+2], h1,h2,h3)

        wf = omega[j]*f1 + omega[j+1]*f2 + omega[j+2]*f3
        dw_dt[j] = f1 - 0.5*cross_x(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j]
        dw_dt[j+1] = f2 - 0.5*cross_y(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+1]
        dw_dt[j+2] = f3 - 0.5*cross_z(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+2]
    end
    return nothing
end

function llg_rhs_gpu(dw_dt::CuArray{T, 1}, m::CuArray{T, 1}, h::CuArray{T, 1},
                 omega::CuArray{T, 1}, alpha::T, gamma::T, precession::Bool, N::Int64) where {T<:AbstractFloat}

    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr llg_rhs_kernal!(dw_dt, m, h, omega, alpha, gamma, precession, N)
    return nothing
end
