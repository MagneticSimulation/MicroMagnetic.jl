function llg_rhs_Cay(dw_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
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

function llg_rhs(dm_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                 alpha::Float64, gamma::Float64, precession::Bool, N::Int64)
  for i = 0:N-1
    j = 3*i+1
    a = -gamma/(1+alpha*alpha)

	f1,f2,f3 = 0.0,0.0,0.0

	mm = m[j]*m[j] + m[j+1]*m[j+1] + m[j+2]*m[j+2]
    mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
    h1 = mm*h[j] - mh*m[j]
    h2 = mm*h[j+1] - mh*m[j+1]
    h3 = mm*h[j+2] - mh*m[j+2]

    if precession
		f1 = cross_x(m[j],m[j+1],m[j+2], h1,h2,h3)
		f2 = cross_y(m[j],m[j+1],m[j+2], h1,h2,h3)
		f3 = cross_z(m[j],m[j+1],m[j+2], h1,h2,h3)
	end

	dm_dt[j] = a * (f1 - h1 * alpha);
    dm_dt[j+1] = a * (f2 - h2 * alpha);
    dm_dt[j+2] = a * (f3 - h3 * alpha);
  end
end

function llg_stt_rhs(dm_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                     h_stt::Array{Float64, 1}, alpha::Float64, beta::Float64, gamma::Float64, N::Int64)
  for i = 0:N-1
    j = 3*i+1
    a = -gamma/(1+alpha*alpha)
    mx, my, mz = m[j], m[j+1], m[j+2]

    mm = mx*mx + my*my + mz*mz
    mh = mx*h[j] + my*h[j+1] + mz*h[j+2]
    h1 = mm*h[j] - mh*mx
    h2 = mm*h[j+1] - mh*my
    h3 = mm*h[j+2] - mh*mz

    f1 = cross_x(mx,my,mz, h1,h2,h3)
    f2 = cross_y(mx,my,mz, h1,h2,h3)
    f3 = cross_z(mx,my,mz, h1,h2,h3)

    dm_dt[j] = a * (f1 - h1 * alpha);
    dm_dt[j+1] = a * (f2 - h2 * alpha);
    dm_dt[j+2] = a * (f3 - h3 * alpha);
    #above is the std LLG equation

    b = 1/(1+alpha*alpha)
    mht = mx*h_stt[j] + my*h_stt[j+1] + mz*h_stt[j+2]
    hp1 = mm*h_stt[j] - mht * mx
    hp2 = mm*h_stt[j+1] - mht * my
    hp3 = mm*h_stt[j+2] - mht * mz

    mth1 = cross_x(mx,my,mz, hp1,hp2,hp3)
    mth2 = cross_y(mx,my,mz, hp1,hp2,hp3)
    mth3 = cross_z(mx,my,mz, hp1,hp2,hp3)

    dm_dt[j] += b * ((1 + alpha* beta) * hp1 - (beta - alpha) * mth1)
    dm_dt[j+1] += b * ((1 + alpha* beta) * hp2 - (beta - alpha) * mth2)
    dm_dt[j+2] += b * ((1 + alpha* beta) * hp3 - (beta - alpha) * mth3)
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

#compute (\vec{u} \cdot \nabla \vec{m}) only for z-direction
function compute_field_stt_trianglar(m::Array{T, 1}, h_stt::Array{T, 1},
                           uz::Array{T, 1}, dz::T,
                           ngbs::Array{Int64, 2}, N::Int64) where {T<:AbstractFloat}

  for i = 1:N
      fx, fy, fz = 0.0, 0.0, 0.0

      #z-direction
      i1 = ngbs[7,i]
      i2 = ngbs[8,i]
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


function llg_stt_rhs_Cay(dw_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                     h_stt::Array{Float64, 1}, omega::Array{Float64, 1}, alpha::Float64,
                     beta::Float64, gamma::Float64, N::Int64)
  for i = 1:N
    j = 3*i-2
    a = gamma/(1+alpha*alpha)
    b = alpha*a
	u = 1.0/(1+alpha*alpha)
    mx, my, mz = m[j], m[j+1], m[j+2]
    hx, hy, hz = h[j], h[j+1], h[j+2]
    hx_stt, hy_stt, hz_stt = h_stt[j], h_stt[j+1], h_stt[j+2]
    ox, oy, oz = omega[j], omega[j+1], omega[j+2]

	mh = mx*hx + my*hy + mz*hz
    hpx = hx - mh*mx
    hpy = hy - mh*my
    hpz = hz - mh*mz

    mht = mx*hx_stt + my*hy_stt + mz*hz_stt
	hpx_stt = hx_stt - mht*mx;
    hpy_stt = hy_stt - mht*my;
    hpz_stt = hz_stt - mht*mz;

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
  llg_rhs_Cay(dw_dt, sim.spin, sim.field, omega, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.nxyz)

  return nothing

end

function llg_call_back(sim::AbstractSim, dm_dt::Array{Float64, 1}, spin::Array{Float64, 1}, t::Float64)

  effective_field(sim, spin, t)
  llg_rhs(dm_dt, spin, sim.field, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.nxyz)

  return nothing
end

function llg_stt_call_back(sim::AbstractSim, dm_dt::Array{Float64, 1}, spin::Array{Float64, 1}, t::Float64)

  driver = sim.driver
  mesh = sim.mesh

  effective_field(sim, spin, t)

  compute_field_stt(spin, driver.h_stt, driver.ux, driver.uy, driver.uz, mesh.dx, mesh.dy, mesh.dz, mesh.ngbs, sim.nxyz)

  llg_stt_rhs(dm_dt, spin, sim.field, driver.h_stt, driver.alpha, driver.beta, driver.gamma, sim.nxyz)

  return nothing

end


function llg_stt_cay_call_back(sim::AbstractSim, dw_dt::Array{Float64, 1}, t::Float64, omega::Array{Float64, 1})

  driver = sim.driver
  mesh = sim.mesh
  omega_to_spin(omega, sim.prespin, sim.spin, sim.nxyz)
  effective_field(sim, sim.spin, t)

  compute_field_stt(sim.spin, driver.h_stt, driver.ux, driver.uy, driver.uz, mesh.dx, mesh.dy, mesh.dz, mesh.ngbs, sim.nxyz)

  llg_stt_rhs_Cay(dw_dt, sim.spin, sim.field, driver.h_stt, omega, driver.alpha, driver.beta, driver.gamma, sim.nxyz)

  return nothing

end
