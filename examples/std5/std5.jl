using JuMag
using NPZ

JuMag.cuda_using_double(true)
function init_fun(i,j,k,dx,dy,dz)
    x = i-10
    y = j-10
    r = (x^2+y^2)^0.5
    if r<2
      return (0,0,1)
    end
    return (y/r, -x/r, 0)
end

mesh_cpu = FDMesh(nx=20,ny=20,nz=2,dx=5e-9,dy=5e-9,dz=5e-9)
mesh_gpu = FDMeshGPU(nx=20,ny=20,nz=2,dx=5e-9,dy=5e-9,dz=5e-9)


function relax_vortex(mesh)
    sim = Sim(mesh,driver="SD",name="std5")

    set_Ms(sim,8e5)
    init_m0(sim,init_fun)
    add_exch(sim,1.3e-11) #exchange length=5.7nm
    add_demag(sim)

    relax(sim,maxsteps=10000,save_m_every=-1)

    if isa(mesh,JuMag.FDMeshGPU)
        name = "m0_gpu"
    else
        name = "m0_cpu"
    end
    save_m(sim,name,npy=true)
end

function applied_current(mesh,ux,beta)

    if isa(mesh,JuMag.FDMeshGPU)
        sim_name = "std5_dyn_gpu"
        m0 = npzread("npys/m0_gpu.npy")
    else
        sim_name = "std5_dyn_cpu"
        m0 = npzread("npys/m0_cpu.npy")
    end
    sim = Sim(mesh,driver="LLG_STT",name=sim_name)
    sim.driver.alpha=0.1
    sim.driver.gamma = 2.211e5
    sim.driver.beta = beta

    set_ux(sim,ux)

    set_Ms(sim,8e5)
    init_m0(sim,m0)
    add_exch(sim,1.3e-11) #exchange length=5.7nm
    add_demag(sim)


    for i=0:800
        run_until(sim,i*1e-11)
        #Rxs,Rys return the x,y coordinate of vortex centre for each layer
        Rxs,Rys = JuMag.compute_guiding_centre(Array(sim.spin),mesh)
        Rx,Ry = sum(Rxs)/2,sum(Rys)/2

        io = open(sim.name*"_centre.txt","a")
        write(io, string(Rx), " ", string(Ry))
        write(io, "\n")
        close(io)

        if i%20 ==0
            println(i)
        end

    end
end

relax_vortex(mesh_cpu)
relax_vortex(mesh_gpu)

#ux=-P*g*mu_B*J/(2*c_e*Ms) where J=1e12 is the current density， P=1 is the polarization rate,
#g is the Lander factor of free electrons, c_e is the charge of an electron.
@info("Start STT!")
applied_current(mesh_cpu,-72.438,0.05)
@info("Start GPU STT!")
applied_current(mesh_gpu,-72.438,0.05)
