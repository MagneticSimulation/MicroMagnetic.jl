function llg_rhs(dw_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                 omega::Array{Float64, 1}, alpha::Float64, gamma::Float64, precession::Bool, N::Int64)
  for i = 0:N-1
    j = 3*i+1
    a = gamma/(1+alpha*alpha)
    b = alpha*a
    mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
    h1 = h[j] - mh*m[j]
    h2 = h[j+1] - mh*m[j+1]
    h3 = h[j+2] - mh*m[j+2]
    f1 = a*h1*precession + b*cross_x(m[j],m[j+1],m[j+2], h1,h2,h3)
    f2 = a*h2*precession + b*cross_y(m[j],m[j+1],m[j+2], h1,h2,h3)
    f3 = a*h3*precession + b*cross_z(m[j],m[j+1],m[j+2], h1,h2,h3)

    wf = omega[j]*f1 + omega[j+1]*f2 + omega[j+2]*f3
    dw_dt[j] = f1 - 0.5*cross_x(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j]
    dw_dt[j+1] = f2 - 0.5*cross_y(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+1]
    dw_dt[j+2] = f3 - 0.5*cross_z(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+2]
  end
end

#compute (\vec{u} \cdot \nabla \vec{m})
function compute_field_stt(m::Array{T, 1}, h_stt::Array{T, 1},
                           ux::Array{T, 1}, uy::Array{T, 1}, uz::Array{T, 1},
                           dx::T, dy::T, dz::T, ngbs::Array{Int64, 2}, N::Int64) where {T<:AbstractFloat}

  for i = 1:N
      fx, fy, fz = 0.0, 0.0, 0.0
      #x-direction
      i1 = ngbs[1,i]
      i2 = ngbs[2,i]
      factor = i1*i2>0 ? 1/(2*dx) : 1/dx
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      j1 = 3*i1-2
      j2 = 3*i2-2
      fx += ux[i] * (m[j2] - m[j1]) * factor;
      fy += ux[i] * (m[j2+1] - m[j1+1]) * factor;
      fz += ux[i] * (m[j2+2] - m[j1+2]) * factor;

      #y-direction
      i1 = ngbs[3,i]
      i2 = ngbs[4,i]
      factor = i1*i2>0 ? 1/(2*dy) : 1/dy
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      j1 = 3*i1-2
      j2 = 3*i2-2
      fx += uy[i] * (m[j2] - m[j1]) * factor;
      fy += uy[i] * (m[j2+1] - m[j1+1]) * factor;
      fz += uy[i] * (m[j2+2] - m[j1+2]) * factor;

      #z-direction
      i1 = ngbs[5,i]
      i2 = ngbs[6,i]
      factor = i1*i2>0 ? 1/(2*dz) : 1/dz
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      j1 = 3*i1-2
      j2 = 3*i2-2
      fx += uz[i] * (m[j2] - m[j1]) * factor;
      fy += uz[i] * (m[j2+1] - m[j1+1]) * factor;
      fz += uz[i] * (m[j2+2] - m[j1+2]) * factor;

      h_stt[3*i-2] = fx
      h_stt[3*i-1] = fy
      h_stt[3*i	] = fz

  end

end


function llg_rhs_stt(dw_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                     h_stt::Array{Float64, 1}, omega::Array{Float64, 1}, alpha::Float64,
					 beta::Float64, gamma::Float64, N::Int64)
  for i = 0:N-1
    j = 3*i+1
    a = gamma/(1+alpha*alpha)
    b = alpha*a
	u = -1.0/(1+alpha*alpha)

	mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
    h1 = h[j] - mh*m[j]
    h2 = h[j+1] - mh*m[j+1]
    h3 = h[j+2] - mh*m[j+2]

    mht = m[j]*h_stt[j] + m[j+1]*h_stt[j+1] + m[j+2]*h_stt[j+2]
	ht1 = h_stt[j] - mht*m[j];
	ht2 = h_stt[j+1] - mht*m[j+1];
	ht3 = h_stt[j+2] - mht*m[j+2];

    f1 = a*h1 + u*(beta-alpha)*ht1
    f2 = a*h2 + u*(beta-alpha)*ht2
    f3 = a*h3 + u*(beta-alpha)*ht3

    ht1 = b*h1 + u*(1+alpha*beta)*ht1
    ht2 = b*h2 + u*(1+alpha*beta)*ht2
    ht3 = b*h2 + u*(1+alpha*beta)*ht3

    f1 += cross_x(m[j],m[j+1],m[j+2], ht1, ht2, ht3)
    f2 += cross_y(m[j],m[j+1],m[j+2], ht1, ht2, ht3)
    f3 += cross_z(m[j],m[j+1],m[j+2], ht1, ht2, ht3)

    wf = omega[j]*f1 + omega[j+1]*f2 + omega[j+2]*f3
    dw_dt[j] = f1 - 0.5*cross_x(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j]
    dw_dt[j+1] = f2 - 0.5*cross_y(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+1]
    dw_dt[j+2] = f3 - 0.5*cross_z(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+2]
  end
end

function llg_call_back(sim::AbstractSim, t::Float64, omega::Array{Float64, 1})

  dw_dt = sim.driver.ode.dw_dt
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)
  llg_rhs(dw_dt, sim.spin, sim.field, omega, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.nxyz)

  return dw_dt

end


function llg_stt_call_back(sim::AbstractSim, t::Float64, omega::Array{Float64, 1})

  driver = sim.driver
  dw_dt = driver.ode.dw_dt
  mesh = sim.mesh
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)
  compute_field_stt(sim.spin, driver.h_stt, driver.ux, driver.uy, driver.uz, mesh.dx, mesh.dy, mesh.dz, mesh.ngbs, sim.nxyz)
  llg_rhs_stt(dw_dt, sim.spin, sim.field, driver.h_stt, omega, driver.alpha, driver.beta, driver.gamma, sim.nxyz)

  return dw_dt

end
