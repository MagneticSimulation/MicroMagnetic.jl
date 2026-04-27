# @title Skyrmion Lattice
# @description Create a skyrmion lattice using parameters from PRL 108 017206 (2012)
# @tags atomistic skyrmion

# Import MicroMagnetic
using MicroMagnetic

# Define Initial Magnetization Function
function m0_fun(i, j, k, dx, dy, dz)
    i0, j0, r = 166, 96, 25
    i1 = i % i0
    j1 = j % j0

    if ((i1 - r)^2 + (j1 - r)^2 < r^2)
        return (0.05, 0.01, -1)
    elseif ((i1 - i0 / 2.0 - r)^2 + (j1 - j0 / 2.0 - r)^2 < r^2)
        return (0.05, 0.01, -1)
    end

    return (0, 0, 1)
end

# Create Mesh
mesh = CubicMesh(; nx=166 * 2, ny=96 * 3, nz=1, pbc="xy")

# Create Simulation
sim = Sim(mesh; driver="SD", name="skx_lattice")

# Set Magnetic Moment
set_mu_s(sim, 1.0)

# Add Interactions
add_exch(sim, 1.0; name="exch")
add_zeeman(sim, (0, 0, 3.75e-3))
add_dmi(sim, 0.09; name="dmi")

# Initialize Magnetization
init_m0(sim, m0_fun)

# Relax System
relax(sim; max_steps=5000, stopping_dmdt=1e-5)