#implement the  Runge-Kutta method with adaptive stepsize using Dormand–Prince pair
#https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Dormand%E2%80%93Prince

#To use this integrator, one needs to provide:
#   - A simulation instance with field "spin" and "prespin"
#   - A call back function
using CuArrays

mutable struct DormandPrinceGPU{T<:AbstractFloat}  <: Integrator
    tol::Float64
    t::Float64
    step::Float64
    step_next::Float64
    facmax::Float64
    facmin::Float64
    safety::Float64
    nsteps::Int64
    nfevals::Int64
    errors::CuArray{T, 1}
    k1::CuArray{T, 1}
    k2::CuArray{T, 1}
    k3::CuArray{T, 1}
    k4::CuArray{T, 1}
    k5::CuArray{T, 1}
    k6::CuArray{T, 1}
    k7::CuArray{T, 1}
    rhs_fun::Function
    succeed::Bool
end

function DormandPrinceGPU(nxyz::Int64, rhs_fun, tol::Float64)
    Float = _cuda_using_double.x ? Float64 : Float32
    errors = CuArrays.zeros(Float,3*nxyz)
    k1 = CuArrays.zeros(Float,3*nxyz)
    k2 = CuArrays.zeros(Float,3*nxyz)
    k3 = CuArrays.zeros(Float,3*nxyz)
    k4 = CuArrays.zeros(Float,3*nxyz)
    k5 = CuArrays.zeros(Float,3*nxyz)
    k6 = CuArrays.zeros(Float,3*nxyz)
    k7 = CuArrays.zeros(Float,3*nxyz)
    facmax = 5.0
    facmin = 0.2
    safety = 0.824
  return DormandPrinceGPU(tol, 0.0, 0.0, 0.0, facmax, facmin, safety, 0, 0, errors,
                k1, k2, k3, k4, k5, k6, k7, rhs_fun, false)
end


#y = y0+c1*k1*step
function dopri5__step1__kernel!(y::CuDeviceArray{T, 1}, y0::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
       @inbounds y[i] = y0[i] + c1*k1[i]*step
    end
	return nothing
end

#y = y0 + (c1*k1+c2*k2)*step
function dopri5__step2__kernel!(y::CuDeviceArray{T, 1}, y0::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                c2::Float64, k2::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = y0[i] + (c1*k1[i]+c2*k2[i])*step
    end
	return nothing
end

function dopri5__step3__kernel!(y::CuDeviceArray{T, 1}, y0::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                c2::Float64, k2::CuDeviceArray{T, 1},
                                c3::Float64, k3::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = y0[i] + (c1*k1[i]+c2*k2[i]+c3*k3[i])*step
    end
	return nothing
end

function dopri5__step4__kernel!(y::CuDeviceArray{T, 1}, y0::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                c2::Float64, k2::CuDeviceArray{T, 1},
                                c3::Float64, k3::CuDeviceArray{T, 1},
                                c4::Float64, k4::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = y0[i] + (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i])*step
    end
	return nothing
end

function dopri5__step5__kernel!(y::CuDeviceArray{T, 1}, y0::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                c2::Float64, k2::CuDeviceArray{T, 1},
                                c3::Float64, k3::CuDeviceArray{T, 1},
                                c4::Float64, k4::CuDeviceArray{T, 1},
                                c5::Float64, k5::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = y0[i] + (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i])*step
    end
	return nothing
end

function dopri5__step6__kernel!(y::CuDeviceArray{T, 1}, y0::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                c2::Float64, k2::CuDeviceArray{T, 1},
                                c3::Float64, k3::CuDeviceArray{T, 1},
                                c4::Float64, k4::CuDeviceArray{T, 1},
                                c5::Float64, k5::CuDeviceArray{T, 1},
                                c6::Float64, k6::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = y0[i] + (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i])*step
    end
	return nothing
end

function dopri5__step7__kernel!(y::CuDeviceArray{T, 1},
                                c1::Float64, k1::CuDeviceArray{T, 1},
                                c2::Float64, k2::CuDeviceArray{T, 1},
                                c3::Float64, k3::CuDeviceArray{T, 1},
                                c4::Float64, k4::CuDeviceArray{T, 1},
                                c5::Float64, k5::CuDeviceArray{T, 1},
                                c6::Float64, k6::CuDeviceArray{T, 1},
                                c7::Float64, k7::CuDeviceArray{T, 1},
                                step::Float64, n::Int64) where {T<:AbstractFloat}
	i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0<i<=n
        @inbounds y[i] = (c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]+c7*k7[i])*step
    end
	return nothing
end

#This function works for GPU+MPI
function dopri5_step_inner_GPU(sim::AbstractSim, step::Float64, t::Float64)
  a = (1/5, 3/10, 4/5, 8/9, 1.0, 1.0)
  b = (1/5, 3/40, 9/40, 44/45, -56/15, 32/9)
  c = (19372/6561, -25360/2187, 64448/6561, -212/729)
  d = (9017/3168, -355/33, 46732/5247, 49/176, -5103/18656)
  v = (35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84)
  w = (71/57600, 0.0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40)
  ode = sim.driver.ode

  N = length(sim.spin)
  blk, thr = CuArrays.cudims(sim.spin)

  sim.spin .= sim.prespin
  ode.rhs_fun(sim, ode.k1, sim.spin, t) #compute k1

  @cuda blocks=blk threads=thr dopri5__step1__kernel!(sim.spin,
                               sim.prespin, b[1], ode.k1, step, N)
  ode.rhs_fun(sim, ode.k2, sim.spin, t + a[1]*step) #k2

  @cuda blocks=blk threads=thr dopri5__step2__kernel!(sim.spin,
                               sim.prespin, b[2], ode.k1,
                               b[3], ode.k2, step, N)
  ode.rhs_fun(sim, ode.k3, sim.spin, t + a[2]*step) #k3

  @cuda blocks=blk threads=thr dopri5__step3__kernel!(sim.spin,
                               sim.prespin,
                               b[4], ode.k1,
                               b[5], ode.k2,
                               b[6], ode.k3, step, N)
  ode.rhs_fun(sim, ode.k4, sim.spin, t + a[3]*step) #k4

  @cuda blocks=blk threads=thr dopri5__step4__kernel!(sim.spin,
                               sim.prespin,
                               c[1], ode.k1,
                               c[2], ode.k2,
                               c[3], ode.k3,
                               c[4], ode.k4, step, N)
  ode.rhs_fun(sim, ode.k5, sim.spin, t + a[4]*step) #k5

  @cuda blocks=blk threads=thr dopri5__step5__kernel!(sim.spin,
                               sim.prespin,
                               d[1], ode.k1,
                               d[2], ode.k2,
                               d[3], ode.k3,
                               d[4], ode.k4,
                               d[5], ode.k5, step, N)
  ode.rhs_fun(sim, ode.k6, sim.spin, t + a[5]*step) #k6

  @cuda blocks=blk threads=thr dopri5__step6__kernel!(sim.spin,
                               sim.prespin,
                               v[1], ode.k1, v[2], ode.k2,
                               v[3], ode.k3, v[4], ode.k4,
                               v[5], ode.k5, v[6], ode.k6,
                               step, N)
  normalise(sim.spin, sim.nxyz) #if we want to copy k7 to k1, we should normalise it here.
  ode.rhs_fun(sim, ode.k7, sim.spin, t + a[6]*step) #k7

  ode.nfevals += 7
  @cuda blocks=blk threads=thr dopri5__step7__kernel!(ode.errors,
                               w[1], ode.k1, w[2], ode.k2,
                               w[3], ode.k3, w[4], ode.k4,
                               w[5], ode.k5, w[6], ode.k6,
                               w[7], ode.k7, step, N)
  abs!(ode.errors)
  max_error =  maximum(ode.errors) + eps()

  return max_error
end


function advance_step(sim::AbstractSim, integrator::DormandPrinceGPU)

    max_nsteps = 100

    t = integrator.t

    sim.prespin .= sim.spin

    if integrator.step_next <= 0
        integrator.step_next = compute_init_step_DP(sim, 1e-12)
        if integrator.step_next<1e-15
            #integrator.step_next = 1e-15
        end
    end

    step_next = integrator.step_next

    nstep = 1
    while true
        max_error = dopri5_step_inner_GPU(sim, step_next, t)/integrator.tol
        if isnan(max_error)
            step_next = 1e-14
        end
        nstep += 1
        if nstep > max_nsteps
            @info("Too many inner steps, the system may not converge!", step_next, max_error)
            return false
        end
        integrator.succeed = (max_error <= 1)

        if integrator.succeed
            integrator.nsteps += 1
            integrator.step = step_next
            integrator.t += integrator.step
            factor =  integrator.safety*(1.0/max_error)^0.2
            integrator.step_next = step_next*min(integrator.facmax, max(integrator.facmin, factor))
            break
        else
            factor =  integrator.safety*(1.0/max_error)^0.25
            step_next = step_next*min(integrator.facmax, max(integrator.facmin, factor))
        end
    end
    return true
end
