# @title Nanobar
# @description Obtain magnetization distribution in a nanobar
# @tags tutorial

# In this example, we consider a nanobar with dimensions 60nm x 10nm x 5 nm
# We used this example to demostrate how to obtain the magnetization distribution using MicroMagnetic.jl

using MicroMagnetic

# Create Mesh
mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);

# Create Simulation
sim = Sim(mesh; driver="SD")

# Set Saturation Magnetization
set_Ms(sim, 8e5)   #Set saturation magnetization Ms=8e5 A/m

# Add Interactions
add_exch(sim, 1e-12);  #Add exchange interaction
add_demag(sim);         #Add demagnetization 

# Initialize Magnetization
init_m0(sim, (1, 1, 0))  #Initialize magnetization

# Relax System
relax(sim; max_steps=2000, stopping_dmdt=0.01)

# Save Magnetization State
save_vtk(sim, "bar"; fields=["exch", "demag"])
