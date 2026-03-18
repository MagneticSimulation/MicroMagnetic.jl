
using MicroMagnetic

using CUDA
using CUDSS

function compute_hysteresis_bs23()
    #mesh = FDMesh(dx=5e-9, dy=5e-9, dz=5e-9, nx=40, ny=40, nz=4)
    mesh = FEMesh("meshes/nanodot.mesh", unit_length=1e-9)

    sim = Sim(mesh, driver="LLG", integrator="BS23", name="disk_fd")
   
    sim.driver.alpha = 0.5
    sim.driver.integrator.tol = 1e-6

    set_Ms(sim, 8e5)

    init_m0(sim, (-1,0,0))
    
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    
    add_zeeman(sim, (0,0,0))
    
    Hs = [i*mT for i=-150:5:150]
    hysteresis(sim, Hs, direction=(1,0,0), full_loop=false, stopping_dmdt=0.05, output="vts")
end

#compute_hysteresis_bs23()

function compute_hysteresis_gpsm()
    
    mesh = FEMesh("meshes/nanodot.mesh", unit_length=1e-9)

    sim = Sim(mesh, driver="LLG", integrator="GPSM", name="disk_fem")
   
    sim.driver.alpha = 0.5
    sim.driver.integrator.step = 1e-12 # 1ps every step

    set_Ms(sim, 8e5)

    init_m0(sim, (-1,0,0))
    
    add_exch(sim, 1.3e-11)
    add_demag(sim)
    
    add_zeeman(sim, (0,0,0))
    
    Hs = [i*mT for i=-150:5:150]
    hysteresis(sim, Hs, direction=(1,0,0), full_loop=false, stopping_dmdt=0.05, output="vts")
end


compute_hysteresis_gpsm()