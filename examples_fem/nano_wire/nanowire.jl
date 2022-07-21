using JuMag
using Test


function init_m_fun(x,y,z)
    if x < 10
        return (-1, 0, 0)
    else
        return (0.8,0.6,0)
    end
end

function relax_nanowire()
    mesh = FEMesh("nanowire.mesh")
    sim = Sim(mesh, driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    add_exch(sim, 1.3e-11)
    add_demag(sim)

    #add_anis(sim, 5.2e5, axis=(0,0.6,0.8))
    relax(sim, maxsteps=2000, stopping_dmdt=0.03, save_vtk_every = -1)
    JuMag.save_inp(sim, "final.inp")

end

relax_nanowire()