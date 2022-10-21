using JuMag
using Test

function  ini_skx(i,j,k,dx,dy,dz)
    i0 = 20
    j0 = 20
    R = 10
    if (i-i0)^2 + (j-j0)^2 <= R^2
        return (0, 0.01, -1)
    else
        return(0, 0.01, 1)
    end
end

function  test_NextExch(Ku, H, vtsname)
    # mesh = SquareMeshGPU(dx=a, dy=a, dz=a, nx=40, ny=40, nz=1, pbc="open")
    mesh = CubicMeshGPU(dx=a, dy=a, dz=a, nx=40, ny=40, nz=1, pbc="open")
    sim = Sim(mesh, driver="SD", name="zss_test/sim")
    # sim.driver.alpha = 0.1
    # sim.driver.gamma = 1

    set_mu_s(sim, 1)
    # init_m0(sim, ini_skx)
    init_m0_random(sim)
    # save_vtk_points(sim, "zss_test/m0.vts")

    add_exch(sim, J1)
    add_next_exch(sim, J2)
    add_next_next_exch(sim, J3)
    add_zeeman(sim, (0,0,H))
    add_anis(sim, Ku)

    relax(sim; maxsteps = 2e4, stopping_dmdt=0.01, save_m_every=5000)
    save_vtk_points(sim, vtsname)

    return
end

if !ispath("zss_test")
    mkpath("zss_test")
end

a = 4e-10
J1 = 1
J2 = -0.8*J1
J3 = -1.2*J1
alpha = 0.1
Ms = 8e5
Ku = 0.1*J1
H = 0.1*J1
vtsname = "zss_test/Ku=1_H=1.vts"
test_NextExch(Ku, H, vtsname)