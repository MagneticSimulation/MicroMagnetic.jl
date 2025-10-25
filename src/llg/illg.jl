"""
GPU kernel to compute the inertial LLG equation,
    dm/dt = v
    dv/dt = -(1/eta) m x v - (gamma/eta) m x (m x H) + (1/tau) m x (m x v) - v^2 m
where eta=alpha*tau.  We rewrite dv/dt to 
    dv/dt = -(1/eta) m x v + (1/eta) m x (m x F) - v^2 m
where F = alpha*v - gamma*H

To solve the equation, we assume y = [m, v].
"""
@kernel function inertial_llg_rhs_kernel!(dy_dt, @Const(y), @Const(h), @Const(pins), alpha::T,
                                 gamma::T, tau::T, N) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds pin::Bool = pins[I]
    inverse_eta::T = 1 / (alpha * tau)

    if pin
        @inbounds dy_dt[j] = 0
        @inbounds dy_dt[j + 1] = 0
        @inbounds dy_dt[j + 2] = 0

        @inbounds dy_dt[j + N] = 0
        @inbounds dy_dt[j + N + 1] = 0
        @inbounds dy_dt[j + N + 2] = 0

    else
        @inbounds mx::T = y[j]
        @inbounds my::T = y[j + 1]
        @inbounds mz::T = y[j + 2]
        
        @inbounds vx::T = y[j + N]
        @inbounds vy::T = y[j + N + 1]
        @inbounds vz::T = y[j + N + 2]

        # F = alpha*v - gamma*H
        @inbounds fx::T = alpha*vx - gamma* h[j]
        @inbounds fy::T = alpha*vy - gamma* h[j+1]
        @inbounds fz::T = alpha*vz - gamma* h[j+2]

        mm::T = mx * mx + my * my + mz * mz
        mf::T = mx * fx + my * fy + mz * fz

        # m x (m x F)
        h1::T = mf * mx - mm * fx
        h2::T = mf * my - mm * fy
        h3::T = mf * mz - mm * fz

        g1::T = cross_x(mx, my, mz, vx, vy, vz)
        g2::T = cross_y(mx, my, mz, vx, vy, vz)
        g3::T = cross_z(mx, my, mz, vx, vy, vz)

        vv::T = vx^2 + vy^2 + vz^2

        # dm/dt = v
        @inbounds dy_dt[j] = vx
        @inbounds dy_dt[j + 1] = vy
        @inbounds dy_dt[j + 2] = vz

        # dv / dt = -(1/eta) m x v + (1/eta) m x (m x F) - v^2 m
        @inbounds dy_dt[j + N] = inverse_eta*(h1 - g1) - vv*mx
        @inbounds dy_dt[j + N + 1] = inverse_eta*(h2 - g2) - vv*my
        @inbounds dy_dt[j + N + 2] = inverse_eta*(h3 - g3) - vv*mz
    end
end

"""
GPU kernel to compute the inertial LLG equation,
    dm/dt = - m x w
    dw/dt = -(1/eta) m x w + (gamma/eta) m x H + (alpha/tau) m x (m x w)
where eta = alpha*tau.  We rewrite dw/dt to 
    dw/dt = (1/eta) m x F 
where F =  gamma H - w + alpha * m x w

To solve the equation, we assume y = [m, w].
"""
@kernel function inertial_llg_rhs_kernel_b!(dy_dt, @Const(y), @Const(h), @Const(pins), alpha::T,
                                 gamma::T, tau::T, N) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * I - 2

    @inbounds pin::Bool = pins[I]
    inverse_eta::T = 1 / (alpha * tau)

    if pin
        @inbounds dy_dt[j] = 0
        @inbounds dy_dt[j + 1] = 0
        @inbounds dy_dt[j + 2] = 0

        @inbounds dy_dt[j + N] = 0
        @inbounds dy_dt[j + N + 1] = 0
        @inbounds dy_dt[j + N + 2] = 0

    else
        @inbounds mx::T = y[j]
        @inbounds my::T = y[j + 1]
        @inbounds mz::T = y[j + 2]
        
        @inbounds wx::T = y[j + N]
        @inbounds wy::T = y[j + N + 1]
        @inbounds wz::T = y[j + N + 2]

        # m x w
        @inbounds mwx::T = cross_x(mx, my, mz, wx, wy, wz)
        @inbounds mwy::T = cross_y(mx, my, mz, wx, wy, wz)
        @inbounds mwz::T = cross_z(mx, my, mz, wx, wy, wz)

        # F = gamma H - w + alpha * m x w
        @inbounds fx::T = gamma * h[j] - wx + alpha*mwx
        @inbounds fy::T = gamma * h[j+1] - wy + alpha*mwy
        @inbounds fz::T = gamma * h[j+2] - wz + alpha*mwz

        # dm/dt = - m x w
        @inbounds dy_dt[j] = -mwx
        @inbounds dy_dt[j + 1] = -mwy
        @inbounds dy_dt[j + 2] = -mwz

        # dw/dt = (1/eta) m x F 
        @inbounds dy_dt[j + N] = inverse_eta*cross_x(mx, my, mz, fx, fy, fz)
        @inbounds dy_dt[j + N + 1] = inverse_eta*cross_y(mx, my, mz, fx, fy, fz)
        @inbounds dy_dt[j + N + 2] = inverse_eta*cross_z(mx, my, mz, fx, fy, fz)
    end
end

"""
LLG call_back function that will be called by the integrator.
"""
function inertial_llg_call_back(sim::AbstractSim, dy_dt::AbstractArray{T,1},
                       y::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    driver = sim.driver

    spin = view(y, 1:(3*N))
    effective_field(sim, spin, t)

    # TODO: For the test in test_illg.jl, the error of inertial_llg_rhs_kernel_b! is much smaller (~8.4e-15),
    # however, it will take 217757 nfevals. As a comparsion, the inertial_llg_rhs_kernel! requires 8800 nfevals
    # and the error is (~8.6e-12). 
    kernel! = inertial_llg_rhs_kernel!(default_backend[], groupsize[])
    kernel!(dy_dt, y, sim.field, sim.pins, T(driver.alpha), T(driver.gamma), T(driver.tau), 3*N, ndrange=N)

    return nothing
end