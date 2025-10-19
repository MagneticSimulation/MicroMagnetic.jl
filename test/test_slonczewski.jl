using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

dot3(m, p) = m[1]*p[1] + m[2]*p[2] + m[3]*p[3]

# H = beta*(epsilon* m x m_p +  xi*m_p)
function analyical(m; J=1e12, tf=2e-9, Ms=8e5, P=0.7, Lambda=2, p=[0,0,1.0], xi=0.05)
   beta = h_bar*J/(mu_0*c_e*tf*Ms)
   epsilon = P*Lambda^2/(Lambda^2+1 + (Lambda^2-1)*dot3(m, p))
   H = beta .* (epsilon .* MicroMagnetic.cross(m, p) .+ xi .* p)
   return H
end

function test_slonczewski_torque(integrator="DormandPrince")
    # Test mesh
    mesh = FDMesh(; nx=1, ny=1, nz=1, dx=1e-9, dy=1e-9, dz=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.tol = 1e-8

    add_zeeman(sim, (0, 0, 1e5))

    # Add Slonczewski STT torque with basic parameters
    t = add_stt(sim, model=:slonczewski, J=1e12, tf=2e-9, Ms=8e5, p=(0,0,1), P=0.7, xi=0.5, Lambda=3.7)

    init_m0(sim, (0.6, 0.8, 0.0), norm=true)
    
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)


    expected = analyical(Array(sim.spin), J=1e12, tf=2e-9, Ms=8e5, p=[0,0,1], P=0.7, xi=0.5, Lambda=3.7)
    @test isapprox(Array(t.field), expected, rtol=1e-6)
end

@using_gpu()
test_functions("SlonczewskiTorque", 
    test_slonczewski_torque
)