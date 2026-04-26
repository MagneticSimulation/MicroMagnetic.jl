# @title Magnetic hopfion
# @description Create and relax a magnetic hopfion
# @tags hopfion tutorial

# Import MicroMagnetic
using MicroMagnetic
@using_gpu()

# Create Mesh
mesh = CubicMesh(nx=64, ny=64, nz=32, dx=1e-9, dy=1e-9, dz=1e-9)

# Define Hopfion Initialization Function
function init_hopfion(x, y, z)
    r = sqrt(x^2+y^2)
    d = 10e-9
    sinf = 2*r*d/(r^2+d^2)
    cosf = sqrt(1-sinf^2)
    mx = x/r*2*sinf*cosf + y*z/r^2*sinf^2
    my = y/r*2*sinf*cosf - x*z/r^2*sinf^2
    mz = cosf^2-sinf^2 + 2*z^2/r^2*sinf^2
    return (mx, my, mz)
end

# Create Simulation
sim = Sim(mesh; driver="SD", name="hopfion")

# Set Magnetic Moment
set_mu_s(sim, 1)

# Initialize Magnetization
init_m0(sim, init_hopfion)

# Add Exchange Interaction
J1 = 1
add_exch(sim, J1; name="exch", J2=-0.164, J3=0, J4=-0.082)

# Relax System
relax(sim; max_steps=50000, stopping_dmdt=1e-3, using_time_factor=false)