function cross_x(a::Array{Float64}, b::Array{Float64})
  return a[2]*b[3]-a[3]*b[2]
end

function cross_y(a::Array{Float64}, b::Array{Float64})
  return a[3]*b[1]-a[1]*b[3]
end

function cross_z(a::Array{Float64}, b::Array{Float64})
  return a[1]*b[2]-a[2]*b[1]
end

function llg_rhs(dm_dt::Array{Float64}, m::Array{Float64}, h::Array{Float64}, alpha::Float64, gamma::Float64)
  n = size(m)[2]
  coeff = -gamma/(1.0 + alpha*alpha)
  hp = [0.0,0.0,0.0]
  mth = [0.0,0.0,0.0]
  for i=1:n
    mm = dot(m[:,i],m[:,i])
    mh = dot(m[:,i],h[:,i])
    hp[:] = mm*h[:,i]-mh*m[:,i]
    mth[1] = cross_x(m[:,i], hp[:])
    mth[2] = cross_y(m[:,i], hp[:])
    mth[3] = cross_z(m[:,i], hp[:])
    dm_dt[:,i] = coeff*(mth[:] - hp[:]*alpha)
    #The above is llg equation
    c = 6*sqrt(dot(dm_dt[:,i], dm_dt[:,i]))

    dm_dt[:,i] += c*(1-mm)*m[:,i]

  end
end
