# ---
# title: Calculating micromagnetic energies
# author: Weiwei Wang
# date: 2022-10-04
# description: an example to demostrate how to compute micromagnetic energies in JuMag.
# tag: tutorial
# ---

# Import JuMag
using JuMag

# Create a finite difference mesh
mesh =  FDMesh(dx=5e-9, dy=5e-9, dz=5e-9, nx=1, ny=1, nz=1)

# Create a Simulation 
sim = Sim(mesh);

# Set the saturation magnetization Ms
set_Ms(sim, 8.6e5)

# Initilize the magnetization to (1,1,1) direction, note that |m|=1 so the magnetization vector m = (1,1,1)/sqrt(3). 
init_m0(sim, (1,1,1))

# Consider the demagnetization energy
demag = add_demag(sim);

# Consider the zeeman energy
zeeman = add_zeeman(sim, (0,0,1e5))

# Calculate the effective field
JuMag.effective_field(sim, sim.spin, 0.0)

# Calculate the energies
println("Demag Energy: ",sum(demag.energy), " J")
println("Zeeman Energy: ",sum(zeeman.energy), " J")
