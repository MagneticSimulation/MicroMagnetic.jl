using MicroMagnetic
using Test
MicroMagnetic.cuda_using_double(true)

function analytical_field_001(mu_s, mx, my, mz, lambda, E, H, omega, delta, t)
  sin_t = sin(omega*t + delta)
  cos_t = cos(omega*t)
  Ex, Ey, Ez = E*sin_t, E*cos_t, 0
  Hx, Hy, Hz = H*cos_t, -H*sin_t, 0

  hx = -(Ez*mx + Ex*mz)*lambda/mu_s + Hx
  hy = (Ez*my + Ey*mz)*lambda/mu_s + Hy
  hz = (-Ex*mx  + Ey*my)*lambda/mu_s + Hz
  return [hx, hy, hz]
end

function analytical_field_110(mu_s, mx, my, mz, lambda, E, H, omega, delta, t)
  sin_t = sin(omega*t + delta)
  cos_t = cos(omega*t)
  Ex, Ey, Ez = E*sin_t, E*cos_t, 0
  Hx, Hy, Hz = H*cos_t, -H*sin_t, 0

  hx = -(Ey*mx + Ex*my)*lambda/mu_s + Hx
  hy = (-Ex*mx + Ez*mz)*lambda/mu_s + Hy
  hz = (Ez*my  + Ey*mz)*lambda/mu_s + Hz
  return [hx, hy, hz]
end

function analytical_field_111(mu_s, mx, my, mz, lambda, E, H, omega, delta, t)
  sin_t = sin(omega*t + delta)
  cos_t = cos(omega*t)
  Ex, Ey, Ez = E*sin_t, E*cos_t, 0
  Hx, Hy, Hz = H*cos_t, -H*sin_t, 0

  lambda = lambda/sqrt(3)
  hx = -(sqrt(2)*Ey*mx + Ez*mx + sqrt(2)*Ex*my + Ex*mz)*lambda/mu_s + Hx
  hy = -(sqrt(2)*Ex*mx + Ez*my + sqrt(2)*Ey*my + Ey*mz)*lambda/mu_s + Hy
  hz = -(Ex*mx  + Ey*my - 2*Ez*mz)*lambda/mu_s + Hz
  return [hx, hy, hz]
end

function test_laser()
    #Test mesh
    mesh =  CubicMeshGPU(nx=1, ny=1, nz=1)

    sim = Sim(mesh, name="spin")

    mu_s = 1.2
    set_mu_s(sim, mu_s)
    
    add_zeeman(sim, (0, 0, 0.1)) #Tesla

    mx, my, mz = 0.8, 0.6, 0
    init_m0(sim, (mx, my, mz))

    lambda = 1.234
    E = 100.0
    B = 2.8
    omega = 8.4
    delta = pi/3
    t = 9.0
    laser = add_magnetoelectric_laser(sim, lambda, E, B, omega, delta=delta, direction=001)

    MicroMagnetic.effective_field(laser, sim, sim.spin, t)
    af001 = analytical_field_001(mu_s, mx, my, mz, lambda, E, B, omega, delta, t)
    @test abs(maximum(af001 - Array(sim.field))) < 1e-10
    
    laser = add_magnetoelectric_laser(sim, lambda, E, B, omega, delta=delta, direction=110, name="l2")
    MicroMagnetic.effective_field(laser, sim, sim.spin, t)
    af110 = analytical_field_110(mu_s, mx, my, mz, lambda, E, B, omega, delta, t)
    @test abs(maximum(af110 - Array(sim.field))) < 1e-10

    laser = add_magnetoelectric_laser(sim, lambda, E, B, omega, delta=delta, direction=111, name="111")
    MicroMagnetic.effective_field(laser, sim, sim.spin, t)
    af111 = analytical_field_111(mu_s, mx, my, mz, lambda, E, B, omega, delta, t)
    @test abs(maximum(af111 - Array(sim.field))) < 1e-10
    
end


test_laser()
