using JuMag
using Printf

JuMag.cuda_using_double(true)

function relax_sim()
    N = 92 #N = 92, 92*5=460, round(3*sqrt(3)*N) = 478
    mesh =  FDMeshGPU(nx=460, ny=478, nz=1, dx=1e-9, dy=1e-9, dz=2e-9, pbc="xy")

    function m0_fun(i,j,k,dx,dy,dz)
        Lx, Ly, bias = N, sqrt(3)*N, N/3
        x, y, r = i%Lx, j%Ly, 10
        m1 = (0.05, 0.01, -1)
        m2 = (0, 0, 1)
        if (x-bias)^2+(y-bias)^2<r^2
            return m1
        elseif (x-Lx/2-bias)^2+(y-Ly/2-bias)^2<r^2
            return m1
        end
        return m2
    end

    name =  @sprintf("N_%d", N)
    sim = Sim(mesh, driver="SD", name=name)
    set_Ms(sim, 3.84e5)

    add_exch(sim, 8.78e-12)
    add_dmi(sim, 1.58e-3)
    add_zeeman(sim, (0,0,180*mT))
    init_m0(sim, m0_fun)

    relax(sim, maxsteps=50000,  stopping_dmdt=0.5, save_vtk_every = -1)
    save_ovf(sim, "skx.ovf")

end

relax_sim()
