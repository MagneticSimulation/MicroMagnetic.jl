using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_torque_time(integrator="DormandPrince")
    mesh = FDMesh(; nx=10, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.gamma = 2.21e5

    amplitude = 1.27e5
    frequency = 2e9  # 2 GHz

    function torque_fun_td(y, m, t)
        y .= 0
        
        y[1:3:end] .= amplitude * sin(2π * frequency * t)      # x: sinusoidal
        y[2:3:end] .= amplitude * cos(2π * frequency * t)      # y: cosinusoidal  
        y[3:3:end] .= amplitude * (1.0 - exp(-t * 1e9))        # z: exponential decay
    end

    t_field = add_torque(sim, torque_fun_td; name="test_torque")

    init_m0(sim, (0.1, 0.2, 0.3), norm=false)
    
    # Test at different times
    test_times = [0.0, 0.45e-9]  # 0, 0.25 ns, 0.5 ns
    
    for t_val in test_times
        
        MicroMagnetic.effective_field(sim, sim.spin, t_val)
        
        # Calculate expected values: (1/γ) m × H(t)
        gamma = 2.21e5
        mx, my, mz = 0.1, 0.2, 0.3
        
        hx = amplitude * sin(2π * frequency * t_val)
        hy = amplitude * cos(2π * frequency * t_val)
        hz = amplitude * (1.0 - exp(-t_val * 1e9))
        
        expected = zeros(length(sim.spin))
        expected[1:3:end] .= (1/gamma) * (my * hz - mz * hy)
        expected[2:3:end] .= (1/gamma) * (mz * hx - mx * hz)
        expected[3:3:end] .= (1/gamma) * (mx * hy - my * hx)
        
        @test isapprox(Array(t_field.field), expected, rtol=1e-6, atol=1e-12)
    end
    
    return true
end

@using_gpu()
test_functions("Torque",
    test_torque_time
)