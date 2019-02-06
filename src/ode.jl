using Sundials

include("field.jl")
include("llg.jl")

type UserData
  x::Float64
end

tt = [0.0]
tout = 0.0
reltol, abstol = 1e-6,1e-6

function cvode_rhs(t, y, ydot, user_data)
  n = 1
  m = reshape(Sundials.asarray(y),3,n)
  dm_dt = reshape(Sundials.asarray(ydot),3,n)
  sim = unsafe_pointer_to_objref(user_data)::SimData
  effective_field(sim, m)
  llg_rhs(dm_dt, m, sim.field, sim.alpha, 1.0)
  return Int32(0)
end

function solver(sim::SimData, ts::Array{Float64})
  n = sim.nxyz
  yout = zeros(3*n)
  y0 = reshape(sim.spin,3*n)
  sim.spin[1] = 1.0

  cvode_mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
  flag = Sundials.CVodeInit(cvode_mem, cvode_rhs, tout, y0)
  flag = Sundials.CVodeSetUserData(cvode_mem, sim)
  flag = Sundials.CVodeSStolerances(cvode_mem, reltol, abstol)
  flag = Sundials.CVSpgmr(cvode_mem, Int64(Sundials.PREC_NONE), 300)
  mx = [1.0]
  my = [0.0]
  for t in ts
    flag = Sundials.CVode(cvode_mem, t, yout, tt, Sundials.CV_NORMAL)
    println("T = ", t, ", Y = ", yout)
    push!(mx, yout[1])
    push!(my, yout[2])
  end
  return mx, my
end
