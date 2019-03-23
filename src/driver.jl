mutable struct EnergyMinimization <: Driver
  gk::Array{Float64, 1}
  tau::Float64
  max_tau::Float64
  min_tau::Float64
  steps::Int64
  nfevals::Int64
end

mutable struct LLG <: Driver
  m::Array{Float64, 1}
  m_next::Array{Float64, 1}
  nsteps::Int64
  nfevals::Int64
end

function create_driver(driver::String, nxyz::Int64) #TODO: FIX ME
    if driver=="SDM"
        gk = zeros(Float64,3*nxyz)
        return EnergyMinimization(gk, 0.0, 1e-4, 1e-14, 0, 0)
    end
    return nothing
end

function cross_product(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64, y3::Float64)
    return (-x3*y2 + x2*y3, x3*y1 - x1*y3, -x2*y1 + x1*y2)
end

function compute_tau(driver::EnergyMinimization, m_pre::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1}, nxyz::Int64)
  if driver.steps == 0
      for i=0:nxyz-1  #compute gk for step 0
        j = 3*i+1
        fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
        gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
        driver.gk[j] = gx
        driver.gk[j+1] = gy
        driver.gk[j+2] = gz
      end
    driver.tau  = driver.min_tau
	#FIXME: check this tau always decrease the system energy?
  else
	sum1, sum2, sum3 = 0.0,0.0,0.0
    for i=0:nxyz-1
        j = 3*i+1
        sx = m[j] - m_pre[j]
        sy = m[j+1] - m_pre[j+1]
        sz = m[j+2] - m_pre[j+2]
        fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
        gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
        fx = gx - driver.gk[j]
        fy = gy - driver.gk[j+1]
        fz = gz - driver.gk[j+2]
        driver.gk[j] = gx
        driver.gk[j+1] = gy
        driver.gk[j+2] = gz
        sum1 += sx*sx+sy*sy+sz*sz
        sum2 += sx*fx+sy*fy+sz*fz
        sum3 += fx*fx+fy*fy+fz*fz
    end
    tau1 = sum2!=0.0 ? sum1/sum2 : 0.0
    tau2 = sum3!=0.0 ? sum2/sum3 : 0.0

	driver.tau = abs(tau1) <= abs(tau2) ? tau2 : tau1
    if driver.tau > driver.max_tau
        driver.tau = sign(driver.tau)*driver.max_tau
    end

  end
  return nothing
end

function run_step(sim::MicroSim, driver::EnergyMinimization)
  effective_field(sim, sim.spin, 0.0)
  compute_tau(driver, sim.prespin, sim.spin, sim.field, sim.nxyz)

  sim.prespin[:] = sim.spin[:]

  h = sim.field
  m = sim.spin
  gk = driver.gk
  tau = driver.tau
  for i=0:sim.nxyz-1
      j = 3*i+1
      fx,fy,fz = cross_product(m[j],m[j+1],m[j+2],h[j],h[j+1],h[j+2])
      factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
      mx = (1-factor)*m[j] - tau*gk[j]
      my = (1-factor)*m[j+1] - tau*gk[j+1]
      mz = (1-factor)*m[j+2] - tau*gk[j+2]
      m[j] = mx/(1+factor)
      m[j+1] = my/(1+factor)
      m[j+2] = mz/(1+factor)
  end
  driver.steps += 1
  max_length_error = error_length_m(sim.spin, sim.nxyz)
  if max_length_error > 1e-15
    normalise(m, sim.nxyz)
  end
  return  nothing
end
