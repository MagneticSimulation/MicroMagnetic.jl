using JuMag
using Test

omega = Array([0.1,0.2,0.3,0,0,0,0.2,0.3,0.4])
spin = Array([0.6,0.8,0,0.6,0.8,0,0.6,0,0.8])
spin_next = zeros(Float64,9)

JuMag.omega_to_spin(omega, spin, spin_next, 3)

expected = Array([0.33816425120772947, 0.9410628019323674, -0.006763285024154589,
                 0.6, 0.8, 0.0, 0.7836829836829837, 0.1361305361305361, 0.6060606060606062])
@test isapprox(spin_next, expected)
