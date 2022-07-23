using JuMag
using Printf

function relax_hexagonal_loop()
    # We first load the mesh generated using netgen with the neutral format.
    mesh = FEMesh("../../test/fem/two_hex.mesh", unit_length=1e-9)

    # We create a simulation with the SD driver. 
    sim = Sim(mesh, driver="SD", name="hexagonal")

    # We set the saturation magnetization of the system.
    set_Ms(sim, 1e6) 

    # We set the initial state of the system.
    init_m0(sim, (1,0,0))

    # In this system, we consider four energies, i.e., exchange, demag, zeeman and anisotropy
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    add_zeeman(sim, (0,0,0))

    #we extract the noramls automatically according to the maximum area of surfaces.
    normals = JuMag.extract_normal_axes_by_maximum_area(mesh)
    @info "debug:", normals
    add_anis(sim, 5e5, axis=normals)

    for i = 500:-20:-500
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

relax_hexagonal_loop()