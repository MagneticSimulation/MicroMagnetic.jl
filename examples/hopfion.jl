# @title Magnetic hopfion
# @description Create and relax a magnetic hopfion
# @tags hopfion tutorial

# Import MicroMagnetic
using MicroMagnetic
@using_gpu()

# Create Mesh
mesh = CubicMesh(nx=64, ny=64, nz=32, dx=1e-9, dy=1e-9, dz=1e-9, pbc="xyz")

# Create Simulation
sim = Sim(mesh; driver="SD", name="hopfion")

# Set Magnetic Moment
set_mu_s(sim, 1)

# Initialize Magnetization
hopf = hopfion(R=20e-9, p=1, q=1)
init_m0(sim, hopf)

# Add Exchange Interaction
J1 = 1
add_exch(sim, J1; name="exch", J2=-0.164, J3=0, J4=-0.082)

# Relax System
relax(sim; max_steps=50000, stopping_dmdt=1e-5, using_time_factor=false)