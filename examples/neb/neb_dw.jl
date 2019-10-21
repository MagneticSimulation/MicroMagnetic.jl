using JuMag
using Printf
#using NPZ

function creat_sim()
    mesh =  FDMesh(nx=200, ny=1, nz=1, dx=4e-9, dy=4e-9, dz=4e-9)
    sim = Sim(mesh, name="neb", driver="SD")
    set_Ms(sim, 8.6e5)

    init_m0(sim, (0.6, 0, -0.8))

    add_exch(sim, 1.3e-11)
    add_anis(sim, 1e5, axis=(1,0,0), name="Kx")
    add_anis(sim, -2e4, axis=(0,0,1), name="Kz")
    return sim
end

function init_m(i,j,k,dx,dy,dz)
    if i<90
        return (-1, 0, 0)
    elseif i>110
        return (1, 0, 0)
    end
    return (0,1,0)
end

sim = creat_sim()

init_images = [(-1, 0, 0), init_m, (1, 0, 0)]

neb = NEB(sim, init_images, [6, 6]; name="dw", driver="LLG")
neb.spring_constant = 1e6
#println(neb.images)

relax(neb, stopping_dmdt=0.1, save_ovf_every=500)
