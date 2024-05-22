using MicroMagnetic
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
    mesh = SquareMeshGPU(dx=a, dy=a, dz=a, nx=40, ny=40, nz=1, pbc="open")
    # mesh = CubicMeshGPU(dx=a, dy=a, dz=a, nx=40, ny=40, nz=1, pbc="open")
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

    relax(sim; maxsteps = 2e4, stopping_dmdt=1.0, save_m_every=2000)
    save_vtk_points(sim, vtsname)

    return
end

if !ispath("zss_test")
    mkpath("zss_test")
end

a = 1
Ms = 8e5
J1 = 1
J2 = -0.8*J1
J3 = -1.2*J1
alpha = 0.1
Ku = 0.1*J1
N_H = 3*meV/(h_bar*2.211e5*Ms) # 自然单位制
N_T = h_bar/(3*meV)
for i = collect(Int, range(21, step=2, stop=29))
    H = i*J1/N_H
    vtsname = "zss_test/Ku=0.1_H=$i.vts"
    test_NextExch(Ku, H, vtsname)
end