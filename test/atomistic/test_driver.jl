using MicroMagnetic
using Test

function init_dw(i, j, k, dx, dy, dz)
    if i < 90
        return (1, 0.1, 0)
    elseif i < 100
        return (0, 1, 0.1)
    else
        return (-1, 0.1, 0)
    end
end

function test_driver()
    #Test mesh
    mesh = CubicMesh(; nx=200, ny=1, nz=1, dx=1e-9)

    sim = Sim(mesh; name="spin", driver="SD")

    set_mu_s(sim, 1.0)

    add_exch(sim, 1.0)
    add_anis(sim, 0.005; axis=(1, 0, 0))

    init_m0(sim, init_dw)
    relax(sim; stopping_dmdt=1e-2, max_steps=1000, using_time_factor=false)

    set_driver(sim; driver="LLG", integrator="DormandPrince")
    sim.driver.alpha = 0.5
    sim.driver.gamma = 1.0

    run_until(sim, 2.0)

    #save_vtk(sim, "dw")

    @test 1 == 1
end

test_driver()
