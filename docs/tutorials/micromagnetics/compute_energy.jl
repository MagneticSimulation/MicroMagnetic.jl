# ---
# title: Calculating micromagnetic energies
# author: Weiwei Wang
# date: 2024-05-16
# description: an example to demostrate how to compute micromagnetic energies in MicroMagnetic.
# tag: tutorial
# ---

# Import MicroMagnetic
using MicroMagnetic

# Create a finite difference mesh
mesh =  FDMesh(dx=5e-9, dy=5e-9, dz=5e-9, nx=10, ny=10, nz=1);

# Create a Simulation with saturation magnetization Ms=8.6e5 A/m
# Initilize the magnetization to (1,1,1) direction, note that |m|=1 so the magnetization vector m = (1,1,1)/sqrt(3). 
sim = create_sim(mesh, Ms=8.6e5, m0=(1,1,1));

# Consider the demagnetization energy
demag = add_demag(sim);

# Consider the zeeman energy
zeeman = add_zeeman(sim, (0,0,1e5));

# Calculate the effective fields and energies
MicroMagnetic.effective_field(sim, sim.spin)

# Calculate the total energies
println("Demag Energy: ",sum(demag.energy), " J")
println("Zeeman Energy: ",sum(zeeman.energy), " J")
