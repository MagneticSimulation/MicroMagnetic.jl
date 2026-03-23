
using MicroMagnetic

function compute_hysteresis_bs23()
    mesh = FEMesh("meshes/nanodot.mesh", unit_length=1e-9)

    sim = Sim(mesh, driver="LLG", integrator="BS23", name="disk")
   
    sim.driver.alpha = 0.5
    sim.driver.integrator.tol = 1e-6

    set_Ms(sim, 8e5)

    init_m0(sim, (-1,0,0))
    
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    
    add_zeeman(sim, (0,0,0))
    
    Hs = [i*mT for i=-150:5:150]
    hysteresis(sim, Hs, direction=(1,0,0), full_loop=false, stopping_dmdt=0.05, output="vtu")
end

#

function compute_hysteresis_gpsm()
    
    mesh = FEMesh("meshes/nanodot.mesh", unit_length=1e-9)

    sim = Sim(mesh, driver="LLG", integrator="GPSM", name="disk")
   
    sim.driver.alpha = 0.5
    sim.driver.integrator.step = 1e-12 # 1ps every step

    set_Ms(sim, 8e5)

    init_m0(sim, (-1,0,0))
    
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    
    add_zeeman(sim, (0,0,0))
    
    Hs = [i*mT for i=-150:5:150]
    hysteresis(sim, Hs, direction=(1,0,0), full_loop=false, stopping_dmdt=0.05, output="vtu")
end

#compute_hysteresis_bs23()
compute_hysteresis_gpsm()

using CairoMakie
fig = plot_ts("disk_llg.txt", x_key="Hx", ["m_x"], x_unit=1/mT, xlabel="H (mT)", ylabel="m", mirror_loop=true);
save("loop.png", fig)