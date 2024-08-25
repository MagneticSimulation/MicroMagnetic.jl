# ---
# title: Skyrmion dynamics under spin transfer torques
# author: Weiwei Wang
# date: 2023-03-01
# description: an example to demostrate how to use MicroMagnetic for skyrmion dynamics
# tag: micromagnetics; skyrmion; spin transfer torque
# ---

# In this example, we will study the skyrmion dynamics in a 2d film. We will save the skyrmion 
# positions to a text file and and generate a movie at the end of simulation. 

# We import the necessary modules
using MicroMagnetic
using Printf

@using_gpu()

set_verbose_logging()

# The studied system is a 800nm x 300nm x 2nm film with periodic boundary conditions.
mesh =  FDMesh(nx=400, ny=150, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy");

# We create a Sim instance using `create_sim` function, and set basic simulation parameters such as
# 'Ms', 'A', 'D' and 'H'.
sim = create_sim(mesh,  Ms=3.87e5, A=5.2e-12, D=1e-3, H=(0,0,160*mT), name="skx");

# We set the initial magnetization configuration to a skyrmion at position (200nm, 80nm) with radius 20 nm.
# Note that the initialized magnetization is a roughly guessing for the skyrimon. 
init_m0_skyrmion(sim, (200e-9, 80e-9), 2e-8)

# We relax the system to obtain the skyrmion profile.
relax(sim, maxsteps=20000, stopping_dmdt=0.01)

# We save the magnetization to vtk, which can be opened using Paraview for 3D visualization. 
MicroMagnetic.save_vtk(sim, "skx")

# After obataining the skyrmion profile, we then move the skyrmion using spin transfer torques.
# So we use change the driver to "LLG_STT" using the `set_driver` function. Meanwhile,
# we set the relevant parameters such as alpha, beta and u. 
set_driver(sim, driver="LLG_STT", alpha=0.05, beta=0.2, ux=-20)

# During the simulation, we need to extract the skyrmion center, so we write a call back function
# in which the skyrmion positions are computed using the `compute_guiding_center` function and 
# saved to a text file with append mode.
function call_back_fun(sim, t)
    Rx, Ry = compute_guiding_center(sim)
    open("XY.txt", "a") do f
        write(f, @sprintf("%g  %g  %g\n", t, Rx, Ry))
    end
end

# We use the `run_sim` function to run the simulation.
# after that, a jld2 file will be created and we can export it to a movie (.mp4) using the `jld2movie`.
if !isfile("assets/skx.mp4")
    run_sim(sim, steps=100, dt=1e-10, save_m_every = 1, call_back=call_back_fun)
    jld2movie("skx.jld2", output="assets/skx.mp4")
end

# ![](./assets/skx.mp4)

