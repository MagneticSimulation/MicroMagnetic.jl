# @title Standard Problem 4
# @description Standard micromagnetic simulation with LLG dynamics
# @tags std tutorial

# Import MicroMagnetic
using MicroMagnetic
@using_gpu()

# Create Mesh
mesh = FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);

# Create Simulation
sim = Sim(mesh; name="std4", driver="SD")

# Set Saturation Magnetization
set_Ms(sim, 8e5)

# Initialize Magnetization
init_m0(sim, (1, 0.25, 0.1))

# Add Exchange Interaction
A = 1.3e-11; # J/m
add_exch(sim, A)

# Add Demagnetization
add_demag(sim)

# Relax System
relax(sim; stopping_dmdt=0.01)

# Set LLG driver
set_driver(sim; driver="LLG", alpha=0.02, gamma=2.211e5)

# Add Zeeman Interaction
add_zeeman(sim, (-24.6mT, 4.3mT, 0))

# Run Simulation
run_sim(sim; steps=100, dt=1e-11, save_m_every=1)