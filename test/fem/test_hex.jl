using JuMag
using Printf

function test_one_hex()
    # We first load the mesh generated using netgen with the neutral format.
    mesh = FEMesh("one_hex.mesh", unit_length=1e-9)

    # We create a simulation with the SD driver. 
    sim = Sim(mesh, driver="SD", name="hexagonal")

    # We set the saturation magnetization of the system.
    set_Ms(sim, 1e6) 

    # We set the initial state of the system.
    init_m0(sim, (1,1,0))

    add_exch(sim, 1.3e-11)
    add_demag(sim)
    add_zeeman(sim, (0,0,1e6))

    normals = JuMag.extract_normal_axes_by_maximum_area(mesh)
    @info "debug:", normals
    add_anis(sim, 5e5, axis=normals)

    relax(sim, maxsteps=10000, stopping_dmdt=0.5)

    JuMag.save_inp(sim, @sprintf("inps/Hx_%d_kA.inp", 10))

end
test_one_hex()