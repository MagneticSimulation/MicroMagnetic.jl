using JuMag
using Printf

function relax()
    mesh =  FDMeshGPU(nx=30, ny=30, nz=30, dx=2e-9, dy=2e-9, dz=2e-9)
    sim = Sim(mesh, name="test", integrator="DormandPrince")
    set_Ms(sim, 8.0e5)
    sim.driver.alpha = 0.1
    sim.driver.gamma = 2.211e5
    sim.driver.ode.tol = 1e-5

    mT = 0.001 / (4*pi*1e-7)
    add_exch(sim, 1.3e-11)
    #add_demag(sim)
    add_zeeman(sim, (-24.6*mT, 4.3*mT, 0))

    T = 500
    add_thermal_noise(sim, T)

    init_m0(sim, (1,1,1))

    for i=1:100

        run_until(sim, 1e-11*i)
        println("Runing ", 1e-11*i, "  dt=", sim.driver.ode.step)
    end
    save_vtk(sim, "test")
end

relax()
