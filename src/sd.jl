#We implement the steepest descent method given by https://aip.scitation.org/doi/10.1063/1.4862839
#and https://doi.org/10.1063/1.4896360  where the Barzilai-Borwein (BB) rule is used to speedup
#the energy minimization

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
	return nothing
  end

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
    tau1 = sum2!=0.0 ? sum1/sum2 : driver.min_tau
    tau2 = sum3!=0.0 ? sum2/sum3 : driver.min_tau
    #println(@sprintf("tau1=%g  tau2=%g", tau1, tau2))

	driver.tau = driver.steps%2 == 0 ? abs(tau2) : abs(tau1)
    if driver.tau > driver.max_tau
        driver.tau = driver.max_tau
    end

  return nothing
end

function run_step(sim::AbstractSim, driver::EnergyMinimization)
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
