#We implement the steepest descent method given by https://aip.scitation.org/doi/10.1063/1.4862839
#and https://doi.org/10.1063/1.4896360  where the Barzilai-Borwein (BB) rule is used to speedup
#the energy minimization
function update_m(m::Array{Float64, 1},h::Array{Float64, 1},tau::Float64,nxyz::Int64)
  n= zeros(3*nxyz)
for i=0:nxyz-1
  j = 3*i+1
  fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
  gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
  factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
  mx = (1-factor)*m[j] - tau*gx
  my = (1-factor)*m[j+1] - tau*gy
  mz = (1-factor)*m[j+2] - tau*gz
  n[j] = mx/(1+factor)
  n[j+1] = my/(1+factor)
  n[j+2] = mz/(1+factor)
end
max_length_error = error_length_m(n, nxyz)
  if max_length_error > 1e-15
    normalise(n, nxyz)
  end
return n
end

function search_tau(sim::AbstractSim)
  driver = sim.driver
  tau1  = driver.min_tau
  tau = 0.0
  tau2 = 0.0
  if isa(sim, MicroSimFEM)
    nxyz = sim.n_nodes
  else
    nxyz = sim.nxyz
  end
  m = sim.spin
  h = sim.field
  dtau = 1e-12
  t = 2

  for i=0:nxyz-1  #compute gk for step 0
    j = 3*i+1
    fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
    gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
    driver.gk[j] = gx
    driver.gk[j+1] = gy
    driver.gk[j+2] = gz
  end
  m1 = update_m(m,h,tau1,nxyz)
  compute_system_energy(sim, m1, 0.0)
  E1 = sum(sim.energy)
  for i = 0:100
    tau2 = tau1 + dtau
    m2 = update_m(m,h,tau2,nxyz)
    compute_system_energy(sim, m2, 0.0)
    E2 = sum(sim.energy)
    if E1>E2
      dtau = t*dtau
      tau = tau1
      E1 = E2
      tau1 = tau2
    elseif i==0
      dtau = -dtau
    else
      break
    end
  end
  a = min(tau,tau2)
  b = max(tau,tau2)
  for j=1:10
    tau1 = a + 1/3*(b-a)
    m1 = update_m(m,h,tau1,nxyz)
    compute_system_energy(sim, m1, 0.0)
    E1 = sum(sim.energy)
    tau2 = a + 2/3*(b-a)
    m2 = update_m(m,h,tau2,nxyz)
    compute_system_energy(sim, m2, 0.0)
    E2 = sum(sim.energy)
    if E1<E2
      a = tau1
    else
      b = tau2
    end
    if b-a<1e-14
      break
    end
  end
  driver.tau = (a+b)/2
  if driver.tau < 0
    driver.tau = -driver.tau
  end
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
  N_spins = isa(sim, MicroSimFEM) ? sim.n_nodes : sim.nxyz

  effective_field(sim, sim.spin, 0.0)
  if driver.steps == 0
    search_tau(sim)
  else
    compute_tau(sim.driver, sim.prespin, sim.spin, sim.field, N_spins)
  end

  sim.prespin[:] = sim.spin[:]

  h = sim.field
  m = sim.spin
  gk = driver.gk
  tau = driver.tau
  
  for i=0:N_spins-1
      if sim.pins[i+1]
          continue
      end
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
  max_length_error = error_length_m(sim.spin, N_spins)
  if max_length_error > 1e-15
    normalise(m, N_spins)
  end
  return  nothing
end
