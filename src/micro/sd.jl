#We implement the steepest descent method given by https://aip.scitation.org/doi/10.1063/1.4862839
#and https://doi.org/10.1063/1.4896360  where the Barzilai-Borwein (BB) rule is used to speedup
#the energy minimization
# m is the magnetization array, in the format [mx1, my1, mz1, ...]
# h is the corresponding effective field
# tau is a scale number
# N is the number of spins, and length(m) == 3*N

@kernel function compute_gk_kernel_1!(gk, @Const(m), @Const(h))
    I = @index(Global)
    j = 3 * I - 2
    @inbounds fx, fy, fz = cross_product(m[j], m[j + 1], m[j + 2], h[j], h[j + 1], h[j + 2])
    @inbounds gx, gy, gz = cross_product(m[j], m[j + 1], m[j + 2], fx, fy, fz)
    @inbounds gk[j] = gx
    @inbounds gk[j + 1] = gy
    @inbounds gk[j + 2] = gz
end

@kernel function compute_gk_kernel_2!(gk, ss, sf, ff, @Const(m), @Const(m_pre), @Const(h))
    I = @index(Global)

    j = 3 * I - 2
    @inbounds sx = m[j] - m_pre[j]
    @inbounds sy = m[j + 1] - m_pre[j + 1]
    @inbounds sz = m[j + 2] - m_pre[j + 2]
    @inbounds fx, fy, fz = cross_product(m[j], m[j + 1], m[j + 2], h[j], h[j + 1], h[j + 2])
    @inbounds gx, gy, gz = cross_product(m[j], m[j + 1], m[j + 2], fx, fy, fz)
    @inbounds fx = gx - gk[j]
    @inbounds fy = gy - gk[j + 1]
    @inbounds fz = gz - gk[j + 2]
    @inbounds gk[j] = gx
    @inbounds gk[j + 1] = gy
    @inbounds gk[j + 2] = gz
    @inbounds ss[I] = sx * sx + sy * sy + sz * sz
    @inbounds sf[I] = sx * fx + sy * fy + sz * fz
    @inbounds ff[I] = fx * fx + fy * fy + fz * fz
end

@kernel function run_step_kernel!(@Const(gk), m, @Const(h), @Const(pins),
                                  tau::T) where {T<:AbstractFloat}
    I = @index(Global)
    @inbounds if !pins[I]
        j = 3 * I - 2
        @inbounds fx, fy, fz = cross_product(m[j], m[j + 1], m[j + 2], h[j], h[j + 1],
                                             h[j + 2])
        factor::T = 0.25 * (fx * fx + fy * fy + fz * fz) * tau^2
        @inbounds mx = (1 - factor) * m[j] - tau * gk[j]
        @inbounds my = (1 - factor) * m[j + 1] - tau * gk[j + 1]
        @inbounds mz = (1 - factor) * m[j + 2] - tau * gk[j + 2]
        @inbounds m[j] = mx / (1 + factor)
        @inbounds m[j + 1] = my / (1 + factor)
        @inbounds m[j + 2] = mz / (1 + factor)
    end
end

function compute_tau(driver::EnergyMinimization, m_pre::AbstractArray{T,1},
                     m::AbstractArray{T,1}, h::AbstractArray{T,1},
                     N::Int64) where {T<:AbstractFloat}
    groupsize = 512

    if driver.steps == 0
        kernel! = compute_gk_kernel_1!(backend[], groupsize)
        kernel!(driver.gk, m, h; ndrange=N)
        KernelAbstractions.synchronize(backend[])

        driver.tau = driver.min_tau
        return nothing
    end

    kernel! = compute_gk_kernel_2!(backend[], groupsize)
    kernel!(driver.gk, driver.ss, driver.sf, driver.ff, m, m_pre, h; ndrange=N)
    KernelAbstractions.synchronize(backend[])

    sum1 = sum(driver.ss) #Is it better to use Float64 for driver.ss?
    sum2 = sum(driver.sf)
    sum3 = sum(driver.ff)

    tau1 = sum2 != T(0) ? sum1 / sum2 : driver.min_tau
    tau2 = sum3 != T(0) ? sum2 / sum3 : driver.min_tau

    driver.tau = driver.steps % 2 == 0 ? abs(tau2) : abs(tau1)
    if driver.tau > driver.max_tau
        driver.tau = driver.max_tau
    end

    return nothing
end

function run_step(sim::AbstractSim, driver::EnergyMinimization)
    effective_field(sim, sim.spin, 0.0)
    compute_tau(driver, sim.prespin, sim.spin, sim.field, sim.n_total)

    sim.prespin .= sim.spin

    groupsize = 512
    kernel! = run_step_kernel!(backend[], groupsize)
    kernel!(driver.gk, sim.spin, sim.field, sim.pins, driver.tau; ndrange=sim.n_total)
    KernelAbstractions.synchronize(backend[])
    driver.steps += 1
    #max_length_error = error_length_m(sim.spin, sim.n_total)
    if driver.steps % 10 == 0
        normalise(sim.spin, sim.n_total)
    end
    return nothing
end
