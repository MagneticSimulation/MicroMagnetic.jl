using JuMag
using Test

JuMag.cuda_using_double(true)

function m0_fun(i, k, nr, nz)

	theta = 2*pi*(i-1)/nr

    return (cos(theta)+0.01,sin(theta),0)
end

function relax_system()
    J = 1*meV
    D = 0.18*J

    mesh = CylindricalTubeMeshGPU(nz=200, nr=30, R=30e-9, dz=2e-9, pbc="z")

    sim = Sim(mesh, driver="SD", name="skx")
    set_mu_s(sim, mu_s_1)

    add_exch(sim, J, name="exch")
    add_dmi(sim, D, name="dmi")

    #Hz= 1*1.5e-2*J / mu_s_1
    #add_zeeman(sim, (0,0,Hz))
	#add_anis_tube(sim, 0.002*J)

    init_m0(sim, m0_fun)

    relax(sim, maxsteps=1000, stopping_dmdt=1e-5, using_time_factor=false)

    JuMag.save_vtu(sim, "skx")

    #Q = compute_skyrmion_number(Array(sim.spin),mesh)
end


relax_system()
