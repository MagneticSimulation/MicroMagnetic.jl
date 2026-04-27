# @title Standard Problem 1 LLG
# @description Standard micromagnetic problem 1 using LLG driver with BS23 integrator
# @tags std hysteresis

# Standard Problem 1, see Fig1 in JAP. 116, 123908 (2014)

using MicroMagnetic
@using_gpu()

# Create Mesh
mesh = FDMesh(dx=10e-9, dy=10e-9, dz=10e-9, nx=100, ny=200, nz=2)

# Create Simulation
sim = Sim(mesh, driver="LLG", integrator="BS23", name="std1")

# Configure LLG Driver
sim.driver.alpha = 0.5
sim.driver.precession = false

# Set Saturation Magnetization
set_Ms(sim, 8e5)

# Initialize Magnetization
init_m0(sim, (-1, 0, 0))

# Add Interactions
add_exch(sim, 1.3e-11)
add_anis(sim, 5e2, axis=(0, 1, 0))
add_demag(sim)

# Add Zeeman Interaction
add_zeeman(sim, (0, 0, 0))

# Run Hysteresis Loop
Hs = [i * mT for i in -50:1:50]
tilt_angle = 1 / 180 * pi
direction = (cos(tilt_angle), sin(tilt_angle), 0)
hysteresis(sim, Hs, direction=direction, full_loop=false, stopping_dmdt=0.1, output="vts")