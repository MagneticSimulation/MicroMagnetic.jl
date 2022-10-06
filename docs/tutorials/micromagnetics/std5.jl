# ---
# title: Micromagnetic Standard Problem 4
# author: Boyao Lyu
# date: 2022-10-02
# description: an example for the micromagnetic standard problem 5 using JuMag.
# link: <https://www.ctcms.nist.gov/~rdm/std5/spec5.xhtml>
# ---

using JuMag
using NPZ

# Enable long float gpu simulation.
JuMag.cuda_using_double(true)

# The system is a rectangular film of magnetic material with dimensions 100 nm × 100 nm × 10 nm.
mesh_cpu = FDMesh(nx=20,ny=20,nz=2,dx=5e-9,dy=5e-9,dz=5e-9)
mesh_gpu = FDMeshGPU(nx=20,ny=20,nz=2,dx=5e-9,dy=5e-9,dz=5e-9)

# Initialize a vortex roughly.
function init_fun(i,j,k,dx,dy,dz)
    x = i-10
    y = j-10
    r = (x^2+y^2)^0.5
    if r<2
      return (0,0,1)
    end
    return (y/r, -x/r, 0)
end

# The material parameters of permalloy
const A = 1.3e-11
const Ms = 8e5

# First step: Relax the system into a energy-optimal state by steepest descent driver.
function relax_vortex(mesh)
    sim = Sim(mesh,driver="SD",name="std5")

    set_Ms(sim, Ms)
    add_exch(sim, A) #exchange length=5.7nm
    add_demag(sim)

    init_m0(sim,init_fun)
    relax(sim,maxsteps=10000,save_m_every=-1)

    npzwrite("std5_m0.npy", Array(sim.spin))
end

# The second step: Run the current driven dynamics.
function applied_current(mesh, ux, beta)
    m0 = npzread("std5_m0.npy")
    sim = Sim(mesh,driver="LLG_STT", name="std5_dyn")
    sim.driver.alpha=0.1
    sim.driver.gamma = 2.211e5
    sim.driver.beta = beta

    set_ux(sim,ux)

    set_Ms(sim,8e5)
    init_m0(sim,m0)
    add_exch(sim,1.3e-11) #exchange length=5.7 nm
    add_demag(sim)

    fname = "std5_center.txt"
    io = open(fname,"w")

    # Use formatstring function to convert the contents into fixed-length string.    
    write(io, JuMag.formatstring("time(s)"))
    write(io, JuMag.formatstring("X(m)"))
    write(io, JuMag.formatstring("Y(m)"))
    close(io)

    for i=0:800
        t = i * 1e-11
        run_until(sim, t)
        # Rxs,Rys return the x,y coordinate of vortex centre for each layer.
        Rxs,Rys = JuMag.compute_guiding_center(Array(sim.spin),mesh)
        Rx,Ry = sum(Rxs)/2,sum(Rys)/2

        io = open(fname,"a")
        write(io, JuMag.formatstring(t))
        write(io, JuMag.formatstring(Rx))
        write(io, JuMag.formatstring(Ry))
        write(io, "\n")
        close(io)

        if i%20 ==0
            println("time=", t)
        end
    end
end

# Choose using either CPU or GPU and run the first step.
mesh = mesh_gpu # mesh = mesh_cpu
isfile("std5_m0.npy") || relax_vortex(mesh)

# Choose a proper current density and run the second step.
# ux=-P*g*mu_B*J/(2*c_e*Ms) 
# where J=1e12 is the current density， P=1 is the polarization rate,
# g is the Lander factor of free electrons, c_e is the charge of an electron.
if !(isfile("std5_dyn_llg_stt.txt") && isfile("std5_center.txt"))
    ux = -72.438
    applied_current(mesh, ux, 0.05)
end

# Plot the initial state.
using CairoMakie
function plot_spatial_m()
    m = npzread("std5_m0.npy")  
    nx, ny, nz = 20, 20, 2    
    m = reshape(m, 3, nx, ny, nz)
  
    fig = Figure(resolution = (1600, 1600))
    ax = Axis(fig[1, 1],
        xlabel = "X (nm)",
        ylabel = "Y (nm)"
        )
    xs = range(0, stop=100, length=nx)
    ys = range(0, stop=100, length=ny)
    mx = m[1, :, :, 1]
    my = m[2, :, :, 1]
    arrows!(ax, xs, ys, mx, my, 
        lengthscale = 2, linewidth = 5, arrowsize = 20, 
        arrowcolor = vec(my), linecolor =  vec(my), align = :center)
  
    save("assets/std5_static.png", fig, px_per_unit = 1) #src  
    return fig
end

plot_spatial_m()

# Now plot the track of vortex centre and the three components.
using DelimitedFiles

function plot_center() 
    data_center = readdlm("std5_center.txt", skipstart=1)
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1],
        xlabel = "dX (m)",
        ylabel = "dY (m)"
    )
    scatterlines!(ax, data_center[:,2] .- 50e-9, data_center[:,3] .- 50e-9, 
        markersize = 3, color = :orange, markercolor = :orange)
    save("assets/std5_center.png", fig) #src
    return fig
end

function plot_m() 
    data = readdlm("std5_dyn_llg_stt.txt", skipstart=2)
    ts, mx, my, mz = data[:,2], data[:,4],data[:,5],data[:,6]
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1],
        xlabel = "t (s)",
        ylabel = "<M>"
    )
    lines!(ax, ts, mx, linewidth=3,label="mx JuMag")
    lines!(ax, ts, my, linewidth=3,label="my JuMag")
    lines!(ax, ts, mz, linewidth=3,label="mz JuMag")

    data_fidimag = readdlm("assets/std5_fidimag_data.txt", skipstart=2)
    ts, mx, my, mz = data_fidimag[:,1], data_fidimag[:,5], data_fidimag[:,6], data_fidimag[:,7]
    lines!(ax, ts, mx, linewidth=3, linestyle=:dash, label="mx fidimag")
    lines!(ax, ts, my, linewidth=3, linestyle=:dash, label="my fidimag")
    lines!(ax, ts, mz, linewidth=3, linestyle=:dash, label="mz fidimag")

    axislegend()
    save("assets/std5_m.png", fig) #src
    return fig
end

plot_center() 
plot_m()