using JuMag
using Printf
using NPZ

function creat_sim()
    mesh =  FDMesh(nx=4, ny=4, nz=4, dx=2.5e-9, dy=2.5e-9, dz=2.5e-9)
    sim = Sim(mesh, name="neb", driver="SD")
    set_Ms(sim, 8e5)

    init_m0(sim, (0.6, 0, -0.8))

    add_exch(sim, 1.3e-11)
    add_anis(sim, 5e4, axis=(1,0,0), name="Kx")
    add_anis(sim, -2e4, axis=(0,0,1), name="Kz")
    return sim
end

function init_m(i,j,k,dx,dy,dz)
    return (1,0,0)
end

creat_sim()

neb = Neb(sim, [ (0,1,1), (-1,0,0)], 6;name="haha")
println(neb.image)

#relax(neb)
