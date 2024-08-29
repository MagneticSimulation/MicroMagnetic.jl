# ---
# title: Standard Problem 5 (sim_with)
# author: Weiwei Wang
# description: an example for the micromagnetic standard problem 5 using MicroMagnetic.
# link: <https://www.ctcms.nist.gov/~rdm/std5/spec5.xhtml>
# ---

# Micromagnetic Standard Problem 5 involves simulating the dynamics of a magnetic vortex in a thin rectangular film under
# the influence of a spin-polarized current. The problem is designed to test the accuracy of micromagnetic simulations in
# capturing the complex behavior of vortex structures, including their motion and deformation. The rectangular film in 
# this example has dimensions of 100 nm × 100 nm × 10 nm, and the simulation tracks the position of the vortex core 
# as it evolves over time due to the combined effects of exchange interactions, demagnetization, and spin-transfer torques.

using MicroMagnetic
using Printf
using CairoMakie

@using_gpu()

# Define a function to initialize a vortex roughly.
function init_fun(i, j, k, dx, dy, dz)
    x = i - 10
    y = j - 10
    r = sqrt(x^2 + y^2)
    if r < 2
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end

# Define a callback function to record vortex center coordinates.
function call_back_fun(sim, t)
    Rx, Ry = compute_guiding_center(sim; z=1)
    open("std5_center.txt", "a") do f
        write(f, @sprintf("%g  %g  %g\n", t, Rx, Ry))
    end
end

# Define simulation parameters.
args = (
    name = "std5",
    task_s = ["relax", "dynamics"],          # List of tasks to perform
    driver_s = ["SD", "LLG_STT"],            # List of drivers to use
    mesh = FDMesh(; nx=20, ny=20, nz=2, dx=5e-9, dy=5e-9, dz=5e-9), # Mesh configuration
    Ms = 8e5,                               # Saturation magnetization
    A = 1.3e-11,                            # Exchange constant
    demag = true,                           # Enable demagnetization
    m0 = init_fun,                          # Initial magnetization function
    alpha = 0.1,                            # Gilbert damping parameter
    ux = -72.438,                           # Effective current density
    beta = 0.05,                            # Nonadiabatic STT parameter
    steps = 160,                            # Number of steps for dynamics
    dt = 0.05ns,                            # Time step size
    stopping_dmdt = 0.01,                   # Stopping criterion for relaxation
    call_back = call_back_fun,              # Callback function for vortex center tracking
    dynamic_m_interval = 1                  # Save magnetization at every step
);

# Run the simulation using `sim_with` function.
sim_with(args);

# Generate a movie for the vortex dynamics.
jld2movie("std5.jld2"; output="assets/std5.mp4", component='z', figsize=(400,400))

# ![](./assets/std5.mp4)
