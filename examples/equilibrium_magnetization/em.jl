using JuMag
using Printf

JuMag.cuda_using_double(true)

function random_m0(i, j, k, nx, ny, nz)
  return 2*rand(3).-1
end

function relax()
    mesh =  FDMeshGPU(nx=64, ny=60, nz=64, dx=2e-9, dy=2e-9, dz=2e-9)
    sim = Sim(mesh, name="test", integrator="DormandPrince")
    set_Ms(sim, 8.0e5)
    sim.driver.alpha = 0.1
    sim.driver.gamma = 2.211e5
    sim.driver.ode.tol = 1e-5

    mT = 0.001 / (4*pi*1e-7)
    add_zeeman(sim, (0, 0, 0))

    add_thermal_noise(sim, 100)

    init_m0(sim, random_m0)

    f = open("M_H.txt", "w")
    write(f, "#H (mT)    <m>   <mz>\n")

    for H = 0:200:2000
        println("Running for $H mT...")
        update_zeeman(sim, (0,0, H*mT))
        for i=1:10000
            JuMag.run_step(sim, sim.driver)
        end
        m = JuMag.average_m(sim)
        ms = sqrt(m[1]^2+m[2]^2+m[3]^2)
        mz = m[3]
        write(f, "$H    $ms   $mz\n")
    end
    close(f)
end

relax()
