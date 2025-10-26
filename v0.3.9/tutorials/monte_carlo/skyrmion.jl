using MicroMagnetic
using NPZ

@using_gpu()  # Enable GPU acceleration if available.

function relax_system(; Hz = 0.1)
    #Create a triangular mesh with periodic boundary conditions in the x and y directions.
    mesh = TriangularMesh(nx = 160, ny = 160, pbc = "xy")

    #Initialize the Monte Carlo simulation object.
    sim = MonteCarlo(mesh; name = "mc")

    #Set up the initial magnetization with random orientation.
    init_m0_random(sim)

    #Add simulation parameters:
    #Exchange interaction.
    add_exch(sim; J = 10.52 * meV)

    #Dzyaloshinskii-Moriya interaction (DMI).
    add_dmi(sim; D = 2.63 * meV, type = "interfacial")

    #Zeeman interaction with external field Hz.
    mu_s = 3.64 * mu_B  # Magnetic moment per spin.
    add_zeeman(sim; Hz = Hz * mu_s)

    #Uniaxial anisotropy.
    add_anis(sim; Ku = 0.29 * meV)

    #Perform high-temperature annealing to prepare the system.
    Ts = [100000, 1000, 500]  # Annealing temperatures (in K).
    for T in Ts
        sim.T = T
        run_sim(sim; max_steps = 10_000, save_vtk_every = -1, save_m_every = -1)
    end

    #Gradual cooling to reach the target temperature of 10K.
    for T in 100:-10:10
        sim.T = T
        run_sim(sim; max_steps = 50_000, save_vtk_every = -1, save_m_every = -1)
    end

    #Save the final results.
    save_vtk(sim, "final.vts")                  # Save magnetization as a VTK file.
    npzwrite("final_m.npy", Array(sim.spin))    # Save magnetization as a NumPy file.
end

#relax_system(Hz = 1.5)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
