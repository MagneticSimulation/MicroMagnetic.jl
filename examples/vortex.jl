# @title Vortex
# @description Create and relax a magnetic vortex
# @tags vortex tutorial

# Create Mesh
mesh = FDMesh(nx=100, ny=100, nz=10, dx=2e-9, dy=2e-9, dz=2e-9);

# Create Simulation
sim = Sim(mesh; name="vortex", driver="SD")

# Set Saturation Magnetization
set_Ms(sim, Cylinder(radius = 100e-9), 8e5)

# Initialize Magnetization
init_m0(sim, vortex(p=1, c=-1))

# Add Exchange Interaction
add_exch(sim, 1.3e-11)

# Add Demagnetization
add_demag(sim)

# Relax System
relax(sim; stopping_dmdt=0.01)