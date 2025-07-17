using MicroMagnetic
using Printf
using CairoMakie

@using_gpu()

function init_fun(i, j, k, dx, dy, dz)
    x = i - 10
    y = j - 10
    r = sqrt(x^2 + y^2)
    if r < 2
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end

args = (name="std5", task_s=["relax", "dynamics"],          # List of tasks to perform
        driver_s=["SD", "LLG_STT"],            # List of drivers to use
        mesh=FDMesh(; nx=20, ny=20, nz=2, dx=5e-9, dy=5e-9, dz=5e-9), # Mesh configuration
        Ms=8e5,                               # Saturation magnetization
        A=1.3e-11,                            # Exchange constant
        demag=true,                           # Enable demagnetization
        m0=init_fun,                          # Initial magnetization function
        alpha=0.1,                            # Gilbert damping parameter
        ux=-72.438,                           # Effective current density
        beta=0.05,                            # Nonadiabatic STT parameter
        steps=160,                            # Number of steps for dynamics
        dt=0.05ns,                            # Time step size
        stopping_dmdt=0.01,                   # Stopping criterion for relaxation
        saver_item=SaverItem(("Rx", "Ry"), ("<m>", "<m>"), compute_guiding_center),    #vortex center tracking
        dynamic_m_interval=1);

sim_with(args);

jld2movie("std5.jld2"; output="assets/std5.mp4", component='z', figsize=(400, 400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
