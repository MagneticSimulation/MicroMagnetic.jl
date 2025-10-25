using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function analytical(t; a=0.3, b=1e18, N=3)
    m = zeros(3*N)
    m[1:3:end] .= cos(a)*sin(b*t^2)
    m[2:3:end] .= sin(a)*sin(b*t^2)
    m[3:3:end] .= cos(b*t^2)
    return m
end

function torque_fun(y, m, t)
    a, b = 0.3, 1e18
    Hx, Hy, Hz = 0, 0, 1e5
    γ, α, tau = 2.21e5, 0.1, 1e-11
    η = α*tau
    sin_bt2 = sin(b * t^2)
    cos_bt2 = cos(b * t^2)
    sin_a = sin(a)
    cos_a = cos(a)
    eta_alpha_t = η + α * t

    Tx = 2b*sin_a*eta_alpha_t + (-γ*Hy + 2b*cos_a*t)*cos_bt2 + γ*Hz*sin_a*sin_bt2
    Ty = (γ*Hx + 2b*sin_a*t)*cos_bt2 - cos_a*(2b*eta_alpha_t + γ*Hz*sin_bt2)
    Tz = (cos_a*γ*Hy - γ*Hx*sin_a - 2b*t)*sin_bt2

    N = length(y)
    z = zeros(eltype(y), N)
    z[1:3:end] .= Tx
    z[2:3:end] .= Ty
    z[3:3:end] .= Tz

    copyto!(y, z)
    return nothing
end

function test_illg(;integrator="DormandPrince")
    mesh = FDMesh(; nx=3, ny=1, nz=1, dx=1e-9)
    sim = Sim(mesh; name="spin", integrator=integrator, driver="InertialLLG")

    set_Ms(sim, 8e5)
    sim.driver.gamma = 2.21e5
    sim.driver.alpha = 0.1
    sim.driver.tau = 1e-11
    sim.driver.integrator.tol = 1e-7
    
    add_zeeman(sim, (0, 0, 1e5))
    add_torque(sim, torque_fun)

    init_m0(sim, (0, 0, 1))
    for i in 1:10
        run_until(sim, 1e-10 * i)
    end

    m_e = analytical(1e-9)
    d = maximum( abs.(m_e .- Array(sim.spin)))
    println(integrator, " max error: ", d)
    error =  eltype(sim.spin) == Float32 ? 2e-6 : 1e-11
    @test d < error
end

@using_gpu()
test_functions("ILLG", test_illg)