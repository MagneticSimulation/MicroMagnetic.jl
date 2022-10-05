#We implement the steepest descent method given by https://aip.scitation.org/doi/10.1063/1.4862839
#and https://doi.org/10.1063/1.4896360  where the Barzilai-Borwein (BB) rule is used to speedup
#the energy minimization

function compute_gk_kernel_1!(gk::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 <index <= nxyz
        j = 3*index - 2
        @inbounds fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
        @inbounds gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
        @inbounds gk[j] = gx
        @inbounds gk[j+1] = gy
        @inbounds gk[j+2] = gz
    end
   return nothing
end

function compute_gk_kernel_2!(gk::CuDeviceArray{T, 1}, ss::CuDeviceArray{T, 1}, sf::CuDeviceArray{T, 1},
                              ff::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, m_pre::CuDeviceArray{T, 1},
                              h::CuDeviceArray{T, 1}, N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0 < index <= N
        j = 3*index - 2
        @inbounds sx = m[j] - m_pre[j]
        @inbounds sy = m[j+1] - m_pre[j+1]
        @inbounds sz = m[j+2] - m_pre[j+2]
        @inbounds fx,fy,fz = cross_product(m[j],m[j+1],m[j+2], h[j],h[j+1],h[j+2])
        @inbounds gx,gy,gz = cross_product(m[j],m[j+1],m[j+2], fx,fy,fz)
        @inbounds fx = gx - gk[j]
        @inbounds fy = gy - gk[j+1]
        @inbounds fz = gz - gk[j+2]
        @inbounds gk[j] = gx
        @inbounds gk[j+1] = gy
        @inbounds gk[j+2] = gz
        @inbounds ss[index] = sx*sx+sy*sy+sz*sz
        @inbounds sf[index] = sx*fx+sy*fy+sz*fz
        @inbounds ff[index] = fx*fx+fy*fy+fz*fz
    end
   return nothing
end

function run_step_kernel!(gk::CuDeviceArray{T, 1}, m::CuDeviceArray{T, 1}, h::CuDeviceArray{T, 1}, pins::CuDeviceArray{Bool, 1}, tau::T, N::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 0<index <= N
        if pins[index]
            return nothing
        end
        j = 3*index - 2
        @inbounds fx,fy,fz = cross_product(m[j],m[j+1],m[j+2],h[j],h[j+1],h[j+2])
        factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
        @inbounds mx = (1-factor)*m[j] - tau*gk[j]
        @inbounds my = (1-factor)*m[j+1] - tau*gk[j+1]
        @inbounds mz = (1-factor)*m[j+2] - tau*gk[j+2]
        @inbounds m[j] = mx/(1+factor)
        @inbounds m[j+1] = my/(1+factor)
        @inbounds m[j+2] = mz/(1+factor)
    end
   return nothing
end


function compute_tau(driver::EnergyMinimizationGPU, m_pre::CuArray{T, 1}, m::CuArray{T, 1},
                     h::CuArray{T, 1}, N::Int64)  where {T<:AbstractFloat}
  if driver.steps == 0
      blk, thr = cudims(N)
      @cuda blocks=blk threads=thr compute_gk_kernel_1!(driver.gk, m, h, N)
      driver.tau  = driver.min_tau
      return nothing
  end

    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr compute_gk_kernel_2!(driver.gk,
                                 driver.ss, driver.sf, driver.ff,
                                 m, m_pre, h, N)
    sum1 = sum(driver.ss) #Is it better to use Float64 for driver.ss?
    sum2 = sum(driver.sf)
    sum3 = sum(driver.ff)

    tau1 = sum2 != T(0) ? sum1/sum2 : driver.min_tau
    tau2 = sum3 != T(0) ? sum2/sum3 : driver.min_tau

	driver.tau = driver.steps%2 == 0 ? abs(tau2) : abs(tau1)
    if driver.tau > driver.max_tau
        driver.tau = driver.max_tau
    end

  return nothing
end

function run_step(sim::AbstractSim, driver::EnergyMinimizationGPU)

  effective_field(sim, sim.spin, 0.0)
  compute_tau(driver, sim.prespin, sim.spin, driver.field, sim.nxyz)

  sim.prespin .= sim.spin

  blk, thr = cudims(sim.nxyz)
  @cuda blocks=blk threads=thr run_step_kernel!(driver.gk, sim.spin,
                               driver.field, sim.pins, driver.tau, sim.nxyz)
  driver.steps += 1
  #max_length_error = error_length_m(sim.spin, sim.nxyz)
  if driver.steps%10 == 0
    normalise(sim.spin, sim.nxyz)
  end
  return  nothing
end
