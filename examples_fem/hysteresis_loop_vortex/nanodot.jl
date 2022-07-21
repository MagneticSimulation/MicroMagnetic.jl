using JuMag
using Printf

#This example simulates a magnetic nanodot, as described in 
#https://nmag.readthedocs.io/en/latest/example_hysteresis_disk/doc.html

function relax_nanodot_loop()
    # We first load the mesh generated using netgen with the neutral format.
    mesh = FEMesh("../nano_disk/nanodisk.mesh", unit_length=1e-9)

    # We create a simulation with the SD driver. 
    sim = Sim(mesh, driver="SD", name="nanodot")

    # We set the saturation magnetization of the system.
    set_Ms(sim, 7.95e5) 

    # We set the initial state of the system.
    init_m0(sim, (1,0,0)) 

    # In this system, we consider three energies, i.e., exchange, demag and zeeman. 
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    add_zeeman(sim, (0,0,0))

    #We traverse different magnetic fields. For each field, we relax the system to obtain 
    #its equilibrium state. 
    for i = 200:-5:-200
        Hx = i*1000.0 # A/m

        #update the zeeman field.
        update_zeeman(sim, (Hx,0,0))

        #The stopping criterion is the stopping_dmdt, typically its value should be in the rangle of [0.01, 1]. 
        relax(sim, maxsteps=10000, stopping_dmdt=0.1)

        #save some simulation data
        save_sim_data(sim)

        #save the magnetization distribution as well
        JuMag.save_inp(sim, @sprintf("inps/Hx_%d_kA.inp", i))
    end

end

relax_nanodot_loop()