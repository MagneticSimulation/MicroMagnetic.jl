#  A Butcher table has the form
#      c1 | a11
#      c2 | a21 a22    
#       . | ...
#      cs | as1 as2 ... ass
#     ----------------------
#         |  b1  b2 ...  bs
#         |  w1  w2 ...  ws

# https://en.wikipedia.org/wiki/Bogacki%E2%80%93Shampine_method
function rk_step_bs23!(sim::AbstractSim, step::Float64, t::Float64, integrator::AdaptiveRK)
    k1, k2, k3, k4 = integrator.ks
    y_current = integrator.y_current
    y_next = integrator.y_next
    y_temp = integrator.y_temp 
    T = eltype(y_current)
    
    # FSAL handling
    if integrator.method.FSAL && integrator.succeed && t >= integrator.t
        k1 .= k4  # k4 is FSAL in BS23
    else
        integrator.rhs_fun(sim, k1, y_current, t)
        integrator.nfevals += 1
    end
    
    # Stage 2: k2 = f(t + 1/2 h, y + 1/2 h k1)
    y_temp .= y_current .+ step * T(0.5) .* k1
    integrator.rhs_fun(sim, k2, y_temp, t + (1/2)*step)
    
    # Stage 3: k3 = f(t + 3/4 h, y + 3/4 h k2)
    y_temp .= y_current .+ step * T(0.75) .* k2
    integrator.rhs_fun(sim, k3, y_temp, t + (3/4)*step)
    
    # Stage 4: k4 = f(t + h, y + 2/9 h k1 + 1/3 h k2 + 4/9 h k3)
    y_next .= y_current .+ step * (T(2/9) .* k1 .+ T(1/3) .* k2 .+ T(4/9) .* k3)
    normalise(y_next, sim.n_total) # Normalization for FSAL stage
    integrator.rhs_fun(sim, k4, y_next, t + step)
    
    integrator.nfevals += 3
        
    y_temp .= -(5/72) .* k1 .+ (1/12) .* k2 .+ (1/9) .* k3 .+ (-1/8) .* k4
    
    max_error = maximum(abs.(y_temp))*step + eps() 
    
    return max_error
end

function rk_step_dopri54!(sim::AbstractSim, step::Float64, t::Float64, I::AdaptiveRK)
    c = (1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0)
    a1 = (1 / 5)
    a2 = (3 / 40, 9 / 40)
    a3 = (44 / 45, -56 / 15, 32 / 9)
    a4 = (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729)
    a5 = (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656)
    a6 = (35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84)
    w = (71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)

    k1, k2, k3, k4, k5, k6, k7 = I.ks
    y_current = I.y_current
    y_temp = I.y_temp  
    y_next = I.y_next
    
    if I.method.FSAL && I.succeed && t >= I.t
        k1 .= k7 # copy k7 directly to k1
    else
        I.rhs_fun(sim, k1, y_current, t) #compute k1
        I.nfevals += 1
    end

    # compute k2
    vector_add2(y_next, y_current, k1, a1[1]*step)
    I.rhs_fun(sim, k2, y_next, t + c[1] * step) 

    # compute k3
    vector_add3(y_next, y_current, k1, k2, a2[1]*step, a2[2]*step)
    I.rhs_fun(sim, k3, y_next, t + c[2] * step) 

    # compute k4
    vector_add4(y_next, y_current, k1, k2, k3, a3[1]*step, a3[2]*step, a3[3]*step)
    I.rhs_fun(sim, k4, y_next, t + c[3] * step) 

    # #compute k5
    vector_add5(y_next, y_current, k1, k2, k3, k4, a4[1]*step, a4[2]*step, a4[3]*step, a4[4]*step)
    I.rhs_fun(sim, k5, y_next, t + c[4] * step) 

    # compute k6
    vector_add6(y_next, y_current, k1, k2, k3, k4, k5, a5[1]*step, a5[2]*step, a5[3]*step, a5[4]*step, a5[5]*step)
    I.rhs_fun(sim, k6, y_next, t + step) 

    # compute k7. note v[2] = 0, so we still use vector_add6
    vector_add6(y_next, y_current, k1, k3, k4, k5, k6, a6[1]*step, a6[3]*step, a6[4]*step, a6[5]*step, a6[6]*step)
    
    normalise(y_next, sim.n_total) #if we want to copy k7 to k1, we should normalise it here.
    
    I.rhs_fun(sim, k7, y_next, t + step) 

    I.nfevals += 6
    vector_add6b(y_temp, k1, k3, k4, k5, k6, k7, w[1], w[3], w[4], w[5], w[6], w[7])

    max_error = maximum(abs.(y_temp))*step + eps()
    
    return max_error
end


# https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
function rk_step_cashkarp54!(sim::AbstractSim, step::Float64, t::Float64, integrator::AdaptiveRK)
    
    # table IV
    a = (0, 1 / 5, 3 / 10, 3 / 5, 1, 7 / 8)
    b2 = (1 / 5)
    b3 = (3 / 40, 9 / 40)
    b4 = (3 / 10, -9 / 10, 6 / 5)
    b5 = (-11 / 54, 5 / 2, -70 / 27, 35 / 27)
    b6 = (1631/55296, 175/512, 575/13824, 44275/110592, 253/4096)
    c = (37/378, 0,	250/621, 125/594, 0, 512/1771)
    w = (-277/64512, 0, 6925/370944, -6925/202752, -277/14336,	277/7084)
    
    I = integrator
    k1, k2, k3, k4, k5, k6 = integrator.ks
    y_current = I.y_current
    y_temp = I.y_temp
    y_next = I.y_next

    #compute k1
    I.rhs_fun(sim, k1, y_current, t) 
    
    # compute k2
    vector_add2(y_next, y_current, k1, b2[1] * step)
    I.rhs_fun(sim, k2, y_next, t + a[2] * step)

    # compute k3
    vector_add3(y_next, y_current, k1, k2, b3[1] * step, b3[2] * step)
    I.rhs_fun(sim, k3, y_next, t + a[3] * step)

    # compute k4
    vector_add4(y_next, y_current, k1, k2, k3, b4[1] * step, b4[2] * step, b4[3] * step)
    I.rhs_fun(sim, k4, y_next, t + a[4] * step)

    # #compute k5
    vector_add5(y_next, y_current, k1, k2, k3, k4, b5[1] * step, b5[2] * step, b5[3] * step,
                b5[4] * step)
    I.rhs_fun(sim, k5, y_next, t + a[5] * step)

    # compute k6
    vector_add6(y_next, y_current, k1, k2, k3, k4, k5, b6[1] * step, b6[2] * step,
                b6[3] * step, b6[4] * step, b6[5] * step)    
    I.rhs_fun(sim, k6, y_next, t + a[6] * step)

    # compute y_next
    vector_add6(y_next, y_current, k1, k3, k4, k5, k6, c[1] * step, c[3] * step, 
                c[4] * step, c[5] * step, c[6] * step)

    normalise(y_next, sim.n_total) 

    I.nfevals += 6
    vector_add5b(y_temp, k1, k3, k4, k5, k6, w[1], w[3], w[4], w[5], w[6])
    max_error = maximum(abs.(y_temp)) * step + eps()

    
    return max_error
end

# https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
function rk_step_fehlberg54!(sim::AbstractSim, step::Float64, t::Float64, integrator::AdaptiveRK)
    # table III
    a = (0, 1 / 4, 3 / 8, 12/13, 1, 1 / 2)
    b2 = (1 / 4)
    b3 = (3/32, 9/32)
    b4 = (1932/2197, -7200/2197, 7296/2197)
    b5 = (439/216, -8,	3680/513, -845/4104)
    b6 = (-8/27, 2,	-3544/2565,	1859/4104, -11/40)
    c = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)
    w = (-1/360, 0, 128/4275, 2197/75240, -1/50, -2/55)

    I = integrator
    k1, k2, k3, k4, k5, k6 = integrator.ks
    y_current = I.y_current
    y_temp = I.y_temp  
    y_next = I.y_next

    #compute k1
    I.rhs_fun(sim, k1, y_current, t) 
    
    # compute k2
    vector_add2(y_next, y_current, k1, b2[1] * step)
    I.rhs_fun(sim, k2, y_next, t + a[2] * step)

    # compute k3
    vector_add3(y_next, y_current, k1, k2, b3[1] * step, b3[2] * step)
    I.rhs_fun(sim, k3, y_next, t + a[3] * step)

    # compute k4
    vector_add4(y_next, y_current, k1, k2, k3, b4[1] * step, b4[2] * step, b4[3] * step)
    I.rhs_fun(sim, k4, y_next, t + a[4] * step)

    # #compute k5
    vector_add5(y_next, y_current, k1, k2, k3, k4, b5[1] * step, b5[2] * step, b5[3] * step,
                b5[4] * step)
    I.rhs_fun(sim, k5, y_next, t + a[5] * step)

    # compute k6
    vector_add6(y_next, y_current, k1, k2, k3, k4, k5, b6[1] * step, b6[2] * step,
                b6[3] * step, b6[4] * step, b6[5] * step)    
    I.rhs_fun(sim, k6, y_next, t + a[6] * step)

    # compute y_next
    vector_add6(y_next, y_current, k1, k3, k4, k5, k6, c[1] * step, c[3] * step, 
                c[4] * step, c[5] * step, c[6] * step)

    normalise(y_next, sim.n_total) 

    I.nfevals += 6
    vector_add5b(y_temp, k1, k3, k4, k5, k6, w[1], w[3], w[4], w[5], w[6])

    max_error = maximum(abs.(y_temp)) * step + eps()
end