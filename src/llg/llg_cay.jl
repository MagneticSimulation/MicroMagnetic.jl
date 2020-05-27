function llg_rhs_Cay(dw_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                 omega::Array{Float64, 1}, pins::Array{Bool, 1},
                 alpha::Float64, gamma::Float64, precession::Bool, N::Int64)
  for i = 0:N-1
      j = 3*i+1
      if pins[i+1]
          dw_dt[j] = 0
          dw_dt[j+1] = 0
          dw_dt[j+2] = 0
          continue
      end

    a = gamma/(1+alpha*alpha)
    b = alpha*a
    mm = m[j]*m[j] + m[j+1]*m[j+1] + m[j+2]*m[j+2]
    mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
    h1 = mm*h[j] - mh*m[j]
    h2 = mm*h[j+1] - mh*m[j+1]
    h3 = mm*h[j+2] - mh*m[j+2]
    f1 = a*h1*precession + b*cross_x(m[j],m[j+1],m[j+2], h1,h2,h3)
    f2 = a*h2*precession + b*cross_y(m[j],m[j+1],m[j+2], h1,h2,h3)
    f3 = a*h3*precession + b*cross_z(m[j],m[j+1],m[j+2], h1,h2,h3)

    wf = omega[j]*f1 + omega[j+1]*f2 + omega[j+2]*f3
    dw_dt[j] = f1 - 0.5*cross_x(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j]
    dw_dt[j+1] = f2 - 0.5*cross_y(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+1]
    dw_dt[j+2] = f3 - 0.5*cross_z(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+2]
  end
end


function llg_stt_rhs_Cay(dw_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                     h_stt::Array{Float64, 1}, omega::Array{Float64, 1}, pins::Array{Bool, 1},
                     alpha::Float64, beta::Float64, gamma::Float64, N::Int64)
  for i = 1:N
    j = 3*i-2
    if pins[i]
        dw_dt[j] = 0
        dw_dt[j+1] = 0
        dw_dt[j+2] = 0
        continue
    end

    a = gamma/(1+alpha*alpha)
    b = alpha*a
	u = 1.0/(1+alpha*alpha)
    mx, my, mz = m[j], m[j+1], m[j+2]
    hx, hy, hz = h[j], h[j+1], h[j+2]
    hx_stt, hy_stt, hz_stt = h_stt[j], h_stt[j+1], h_stt[j+2]
    ox, oy, oz = omega[j], omega[j+1], omega[j+2]

    mm = mx*mx + my*my + mz*mz
	mh = mx*hx + my*hy + mz*hz
    hpx = mm*hx - mh*mx
    hpy = mm*hy - mh*my
    hpz = mm*hz - mh*mz

    mht = mx*hx_stt + my*hy_stt + mz*hz_stt
	hpx_stt = mm*hx_stt - mht*mx;
    hpy_stt = mm*hy_stt - mht*my;
    hpz_stt = mm*hz_stt - mht*mz;

    fx = a*hx + u*(beta-alpha)*hpx_stt
    fy = a*hy + u*(beta-alpha)*hpy_stt
    fz = a*hz + u*(beta-alpha)*hpz_stt

    ht1 = b*hx + u*(1+alpha*beta)*hpx_stt
    ht2 = b*hy + u*(1+alpha*beta)*hpy_stt
    ht3 = b*hz + u*(1+alpha*beta)*hpz_stt

    fx += cross_x(mx, my, mz, ht1, ht2, ht3)
    fy += cross_y(mx, my, mz, ht1, ht2, ht3)
    fz += cross_z(mx, my, mz, ht1, ht2, ht3)

    wf = ox*fx + oy*fy + oz*fz
    dw_dt[j] = fx - 0.5*cross_x(ox, oy, oz, fx, fy, fz) + 0.25*wf*ox
    dw_dt[j+1] = fy - 0.5*cross_y(ox, oy, oz, fx, fy, fz) + 0.25*wf*oy
    dw_dt[j+2] = fz - 0.5*cross_z(ox, oy, oz, fx, fy, fz) + 0.25*wf*oz
  end
end

function llg_cay_call_back(sim::AbstractSim, dw_dt::Array{Float64, 1}, t::Float64, omega::Array{Float64, 1})

  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)
  llg_rhs_Cay(dw_dt, sim.spin, sim.field, omega, sim.pins, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.nxyz)

  return nothing

end


function llg_stt_cay_call_back(sim::AbstractSim, dw_dt::Array{Float64, 1}, t::Float64, omega::Array{Float64, 1})

  driver = sim.driver
  mesh = sim.mesh
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)

  compute_field_stt(sim.spin, driver.h_stt, driver.ux, driver.uy, driver.uz, mesh.dx, mesh.dy, mesh.dz, mesh.ngbs, sim.nxyz)

  llg_stt_rhs_Cay(dw_dt, sim.spin, sim.field, driver.h_stt, omega, sim.pins, driver.alpha, driver.beta, driver.gamma, sim.nxyz)

  return nothing

end
