using JuMag
using Printf

#This example simulates an ellipsoidal magnetic particle. As the particle size is relatively small, 
#the Stoner-Wohlfarth model can very well describe the magnetic particle's behavior under external fields.
#The parameters used in this example are the same as https://nmag.readthedocs.io/en/latest/example_stoner_wohlfarth/doc.html

function relax_sw_loop()
    # We first load the mesh generated using netgen with the neutral format.
    mesh = FEMesh("ellipsoid.mesh")

    # We create a simulation with the SD driver. 
    sim = Sim(mesh, driver="SD", name="sw")

    # We set the saturation magnetization of the system.
    set_Ms(sim, 1e6) 

    # We set the initial state of the system.
    init_m0(sim, (-1,-1,0)) 

    # In this system, we consider three energies, i.e., exchange, demag and zeeman. 
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    add_zeeman(sim, (0,0,0))

    #We traverse different magnetic fields. For each field, we relax the system to obtain 
    #its equilibrium state. 
    for i = 500:-10:-500
        H = i*1000.0 # A/m
        Hx = H*sqrt(2)/2
        Hy = H*sqrt(2)/2

        #update the zeeman field.
        update_zeeman(sim, (Hx,Hy,0))

        #The stopping criterion is the stopping_dmdt, typically its value should be in the rangle of [0.01, 1]. 
        relax(sim, maxsteps=10000, stopping_dmdt=0.05)

        #save some simulation data
        save_sim_data(sim)

        #save the magnetization distribution as well
        JuMag.save_inp(sim, @sprintf("inps/H_%d_kA.inp", i))
    end

end

relax_sw_loop()