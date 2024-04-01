# We implement the steepest descent method given by https://aip.scitation.org/doi/10.1063/1.4862839
# and https://doi.org/10.1063/1.4896360  where the Barzilai-Borwein (BB) rule is used to speedup
# the energy minimization
# m is the magnetization array, in the format [mx1, my1, mz1, ...]
# h is the corresponding effective field
# tau is a scale number
# N is the number of spins, and length(m) == 3*N
function update_m(m::Array{Float64, 1}, h::Array{Float64, 1}, tau::Float64, N::Int64)
  n = zeros(3*N)
  for i=0:N-1
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
  max_length_error = error_length_m(n, N)
  if max_length_error > 1e-15
    normalise(n, N)
  end
return n
end

function search_tau(sim::AbstractSim)
  driver = sim.driver
  tau1  = driver.min_tau
  tau = 0.0
  tau2 = 0.0

  N = sim.n_total
  
  m = sim.spin
  h = sim.field
  dtau = 1e-12
  t = 2

  for i=0:N-1  #compute gk for step 0
    j = 3*i+1
    fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
    gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
    driver.gk[j] = gx
    driver.gk[j+1] = gy
    driver.gk[j+2] = gz
  end
  m1 = update_m(m,h,tau1,N)
  compute_system_energy(sim, m1, 0.0)
  E1 = sum(sim.energy)
  for i = 0:100
    tau2 = tau1 + dtau
    m2 = update_m(m,h,tau2,N)
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
    m1 = update_m(m,h,tau1,N)
    compute_system_energy(sim, m1, 0.0)
    E1 = sum(sim.energy)
    tau2 = a + 2/3*(b-a)
    m2 = update_m(m,h,tau2,N)
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



function compute_tau(driver::EnergyMinimization, m_pre::Array{Float64, 1}, m::Array{Float64, 1}, h::Array{Float64, 1}, N::Int64)
  if driver.steps == 0
      for i=0:N-1  #compute gk for step 0
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

    sum1 = zeros(Threads.nthreads())
    sum2 = zeros(Threads.nthreads())
    sum3 = zeros(Threads.nthreads())
    #sum1, sum2, sum3 = 0.0,0.0,0.0
    Threads.@threads for i=0:N-1
        local j = 3*i+1
        local sx = m[j] - m_pre[j]
        local sy = m[j+1] - m_pre[j+1]
        local sz = m[j+2] - m_pre[j+2]
        local fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
        local gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
        fx = gx - driver.gk[j]
        fy = gy - driver.gk[j+1]
        fz = gz - driver.gk[j+2]
        driver.gk[j] = gx
        driver.gk[j+1] = gy
        driver.gk[j+2] = gz
        sum1[Threads.threadid()] += sx*sx+sy*sy+sz*sz
        sum2[Threads.threadid()] += sx*fx+sy*fy+sz*fz
        sum3[Threads.threadid()] += fx*fx+fy*fy+fz*fz
    end

    sum1, sum2, sum3 = sum(sum1), sum(sum2), sum(sum3)
    #@info @sprintf("sum1=%g  sum2=%g  sum3=%g",sum1, sum2, sum3)
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

  N_spins = sim.n_total

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
  
  Threads.@threads for i=0:N_spins-1
      if sim.pins[i+1]
          continue
      end
      local j = 3*i+1
      local fx,fy,fz = cross_product(m[j],m[j+1],m[j+2],h[j],h[j+1],h[j+2])
      local factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
      local mx = (1-factor)*m[j] - tau*gk[j]
      local my = (1-factor)*m[j+1] - tau*gk[j+1]
      local mz = (1-factor)*m[j+2] - tau*gk[j+2]
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
