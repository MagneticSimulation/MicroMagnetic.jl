using JuMag
using Test

function random_fun(x, y, z)
    return 2*rand(3) .- 1
end

function relax_nanodisk()
    mesh = FEMesh("nanodisk.mesh")
    sim = Sim(mesh, driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, random_fun)

    add_exch(sim, 1.3e-11)
    add_demag(sim)

    #add_anis(sim, 5.2e5, axis=(0,0.6,0.8))
    JuMag.save_inp(sim, "init.inp")
    relax(sim, maxsteps=2000, stopping_dmdt=0.05, save_vtk_every = -1)
    JuMag.save_inp(sim, "final.inp")

end

relax_nanodisk()