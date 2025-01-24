
using MicroMagnetic

# Enable GPU acceleration
@using_gpu()

# Define the system geometry: a film with thickness t = 3 nm, length L = 500 nm, and width d = 125 nm.
# Gather all the parameters related to standard problem 4:
args = (
    name = "std4",
    task_s = ["relax", "dynamics"],           # List of tasks
    mesh = FDMesh(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9),
    Ms = 8e5,                                 # Saturation magnetization
    A = 1.3e-11,                              # Exchange constant
    demag = true,                             # Enable demagnetization
    m0 = (1, 0.25, 0.1),                      # Initial magnetization
    alpha = 0.02,                             # Gilbert damping
    steps = 100,                              # Number of steps for dynamics
    dt = 0.01ns,                              # Time step size
    stopping_dmdt = 0.01,                     # Stopping criterion for relaxation
    H_s = [(0,0,0), (-24.6mT, 4.3mT, 0)]      # Static field sweep
);

sim = sim_with(args);
