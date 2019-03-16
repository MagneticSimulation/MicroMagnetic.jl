function cross_x(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64, y3::Float64)
    return -x3*y2 + x2*y3
end

function cross_y(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64, y3::Float64)
    return x3*y1 - x1*y3
end

function cross_z(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64, y3::Float64)
    return -x2*y1 + x1*y2
end

function llg_rhs(dw_dt::Array{Float64}, m::Array{Float64}, h::Array{Float64},
                 omega::Array{Float64}, alpha::Float64, gamma::Float64, precession::Bool, N::Int64)
  for i = 0:N-1
    j = 3*i+1
    a = -gamma/(1+alpha*alpha)
    b = alpha*a
    mh = m[j]*h[j] + m[j+1]*h[j+1] + m[j+2]*h[j+2]
    h1 = h[j] - mh*m[j]
    h2 = h[j+1] - mh*m[j+1]
    h3 = h[j+2] - mh*m[j+2]
    f1 = -a*h1*precession - b*cross_x(m[j],m[j+1],m[j+2], h1,h2,h3)
    f2 = -a*h2*precession - b*cross_y(m[j],m[j+1],m[j+2], h1,h2,h3)
    f3 = -a*h3*precession - b*cross_z(m[j],m[j+1],m[j+2], h1,h2,h3)

    wf = omega[j]*f1 + omega[j+1]*f2 + omega[j+2]*f3
    dw_dt[j] = f1 - 0.5*cross_x(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j]
    dw_dt[j+1] = f2 - 0.5*cross_y(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+1]
    dw_dt[j+2] = f3 - 0.5*cross_z(omega[j], omega[j+1], omega[j+2], f1, f2, f3) + 0.25*wf*omega[j+2]
  end
end
