"""
Compute the LLG equation without precession term
    dm/dt =  - alpha*gamma_L* m x (m x H)
where gamma_L = gamma/(1+alpha^2).
"""
function llg_rhs(dm_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1},
                 alpha::Float64, gamma::Float64, N::Int64)
  for i = 0:N-1
    j = 3*i+1

    a = -gamma/(1+alpha*alpha)

	f1,f2,f3 = 0.0,0.0,0.0

	mm = m[j]*m[j] + m[j+1]*m[j+1] + m[j+2]*m[j+2]
    mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
    h1 = mm*h[j] - mh*m[j]
    h2 = mm*h[j+1] - mh*m[j+1]
    h3 = mm*h[j+2] - mh*m[j+2]

	dm_dt[j] = a * (f1 - h1 * alpha);
    dm_dt[j+1] = a * (f2 - h2 * alpha);
    dm_dt[j+2] = a * (f3 - h3 * alpha);
  end
end

"""
Compute the standard LLG equation,
    dm/dt = - gamma_L * (m x H) - alpha*gamma_L* m x (m x H)
where gamma_L = gamma/(1+alpha^2).
"""
function llg_rhs(dm_dt::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1}, pins::Array{Bool, 1},
                 alpha::Float64, gamma::Float64, precession::Bool, N::Int64)
  for i = 0:N-1
    j = 3*i+1
    if (pins[i+1])
        dm_dt[j] = 0;
        dm_dt[j+1] = 0;
        dm_dt[j+2] = 0;
        continue;
    end

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

"""
Compute STT torque and add it to dm/dt:
    dm/dt += (1+alpha*beta)/(1+alpha^2)*tau - (beta-alpha)/(1+alpha^2)*(m x tau)
where tau = (u.nabla) m, here tau is represented by h_stt.
"""
function add_stt_rhs(dm_dt::Array{Float64, 1}, m::Array{Float64, 1}, h_stt::Array{Float64, 1}, pins::Array{Bool, 1},
                    alpha::Float64, beta::Float64, N::Int64)

  c1 = (1 + alpha* beta)/(1+alpha*alpha)
  c2 = (beta-alpha)/(1+alpha*alpha)

  for i = 0:N-1
    j = 3*i+1
    if (pins[i+1])
        dm_dt[j] = 0;
        dm_dt[j+1] = 0;
        dm_dt[j+2] = 0;
        continue;
    end

    mx, my, mz = m[j], m[j+1], m[j+2]
    mm = mx*mx + my*my + mz*mz
    mht = mx*h_stt[j] + my*h_stt[j+1] + mz*h_stt[j+2]

    #H_perp = (m.m)H - (m.H)m
    hp1 = mm*h_stt[j] - mht * mx
    hp2 = mm*h_stt[j+1] - mht * my
    hp3 = mm*h_stt[j+2] - mht * mz

    mth1 = cross_x(mx,my,mz, hp1,hp2,hp3)
    mth2 = cross_y(mx,my,mz, hp1,hp2,hp3)
    mth3 = cross_z(mx,my,mz, hp1,hp2,hp3)

    dm_dt[j] += c1 * hp1 - c2 * mth1
    dm_dt[j+1] += c1 * hp2 - c2 * mth2
    dm_dt[j+2] += c1 * hp3 - c2 * mth3
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

"""
LLG call_back function that will be called by the integrator.
"""
function llg_call_back(sim::AbstractSim, dm_dt::Array{Float64, 1}, spin::Array{Float64, 1}, t::Float64)

  effective_field(sim, spin, t)
  llg_rhs(dm_dt, spin, sim.field, sim.pins, sim.driver.alpha, sim.driver.gamma, sim.driver.precession, sim.n_nodes)

  return nothing
end

"""
LLG with STT extrension call_back function that will be called by the integrator.
"""
function llg_stt_call_back(sim::AbstractSim, dm_dt::Array{Float64, 1}, spin::Array{Float64, 1}, t::Float64)

  driver = sim.driver
  mesh = sim.mesh

  effective_field(sim, spin, t)
  compute_field_stt(spin, driver.h_stt, driver.ux, driver.uy, driver.uz, mesh.dx, mesh.dy, mesh.dz, mesh.ngbs, sim.n_nodes)

  llg_rhs(dm_dt, spin, sim.field, sim.pins, sim.driver.alpha, sim.driver.gamma, true, sim.n_nodes)

  add_stt_rhs(dm_dt, spin, driver.h_stt, sim.pins, driver.alpha, driver.beta, sim.n_nodes)

  return nothing

end
