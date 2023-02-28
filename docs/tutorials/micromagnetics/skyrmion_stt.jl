using JuMag
JuMag.cuda_using_double(true)

function init_skx(i,j,k,dx,dy,dz)
    if (i-100)^2 + (j-80)^2 < 10^2
        return (0,0.01,-1)
    end
    return (0, 0, 1)
end

mesh =  FDMeshGPU(nx=300, ny=100, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")

sim = create_sim(mesh, A = 5.2e-12, D = 1e-3, H=(0,0,160*mT), Ms=3.87e5, m0=init_skx, name="skx")

relax(sim, maxsteps=20000,  stopping_dmdt=0.01)

JuMag.save_vtk_points(sim, "skx")

set_driver(sim, driver="LLG_STT", alpha=0.2, beta=0.05, ux=1)

run_sim(sim, steps=20, dt=1e-11, save_m_every = 1)

#save_vtk(sim, "mt")

jdl2avi("skx.jdl2")


