# ---
# title: Skyrmion motion under current
# author: Weiwei Wang
# date: 2022-10-05
# description: to repeat the figure 2 in Nat Commun **4** 1463 (2013).
# tag: atomistic; skyrmion
# ---

#=====================
In JuMag, the current density is denoted using u
```math
\mathbf{u} = \frac{p a^3}{2 e S} \mathbf{j}.
```

!!! note "Used parameters in the simulation"
    |Parameter | Value  | 
    | :----:   | :----: | 
    | Lattice constant | $a = 0.5$ nm |
    | Spin length      | $S = 1$      | 
    | Magnetic moment  |  $\mu_s = 2 \mu_B$ |
    | Excahnge constant |  $J = 1$ meV   |
    | DMI         | $D/J = 0.18$  |     |
    | External field  | $H \mu_s /J  = 0.015$ | 
    | Spin-polarization rate        |  $p = 0.2$ |
    | Conversion factor u2j |  $1.282 \times 10^{10}$ |
    | Gyromagnetic ratio    |  $\gamma=2 \mu_B/\hbar$
=====================#

using JuMag
using Printf
using NPZ
using DelimitedFiles
using CairoMakie

# We define a function to roughly initialize a skyrmion
function m0_fun(i, j, k, dx, dy, dz)
    r = 25
    if ((i - 80)^2 + (j - 100)^2 < r^2)
        return (0.05, 0.01, -1)
    end
    return (0,0,1)
end

# We create a function to describe the basic setup of the problem
function basic_setup(;driver="SD", m0=(0,0,1))

    mesh = CubicMesh(nx=500, ny=200, nz=1, dx=0.5e-9, dy=0.5e-9, pbc="xy")

    sim = Sim(mesh, driver=driver, name="skx")
    set_mu_s(sim, mu_s_1)

    init_m0(sim, m0)
    
    J = 1*meV
    add_exch(sim, J, name="exch")

    D = 0.18*J
    add_dmi(sim, D, name="dmi")
    
    Hz= 1.5e-2*J / mu_s_1
    add_zeeman(sim, (0,0,Hz))

    return sim

end 

# We define a function to relax the system to get a skyrmion
function relax_system()
    sim = basic_setup(m0=m0_fun)
    
    relax(sim, maxsteps=5000, stopping_dmdt=0.001, using_time_factor=false)

    npzwrite("skx.npy", Array(sim.spin))

    save_vtk(sim, "skx")
end

# We relax the system to see whether the skyrmion is created.
if !isfile("skx.npy")
    relax_system()
end


# We define a function to move the skyrmion using stt
function drive_skyrmion(;u=2, alpha=0.1, beta=0.2, dt=1e-10, steps=10)

    sim = basic_setup(m0=npzread("skx.npy"), driver="LLG_STT")
    sim.driver.alpha = alpha
    sim.driver.beta = beta
    sim.driver.gamma = 2*mu_B/h_bar

    set_ux(sim, -u)
    
    isdir("assets") || mkdir("assets")
    f = open(@sprintf("assets/pos_u_%g_alpha_%g_beta_%g.txt", u, alpha, beta), "w")
    for i = 0:steps
        t = dt*i
        println("running at $t")
        run_until(sim, t, save_data=false)
        Rx, Ry = compute_guiding_center(sim)
        write(f, "$t $Rx $Ry\n")
    end
    close(f)

    #save_vtk(sim, @sprintf("skx_u_%g", u))
end

function plot_XY(filename)

    data = readdlm(filename)
    ts, X, Y = data[:,1]*1e9, data[:,2]*1e9,  data[:,3]*1e9
    
    fig = Figure(resolution = (800, 480), fontsize=28)
    ax = Axis(fig[1, 1],
        xlabel = "Time (ns)",
        ylabel = "Displacement (nm)"
    )

    lines!(ax, ts, X.-X[1], label="X")
    lines!(ax, ts, Y.-Y[1], label="Y")
    axislegend()

    save("assets/xy.png",fig)
    return fig

end

# We run the simulation to see how the skyrmion moves, and 
# we find that its trajectory is a straight line
u=2; alpha=0.1; beta=0.2
name = @sprintf("assets/pos_u_%g_alpha_%g_beta_%g.txt", u, alpha, beta)
if !isfile(name)
    drive_skyrmion(u=u, alpha=alpha, beta=beta)
end
plot_XY(name)

# run a series of u
us = [i for i = 1:6]
for u in us
    local name = @sprintf("assets/pos_u_%g_alpha_%g_beta_%g.txt", u, alpha, beta)
    isfile(name) || drive_skyrmion(u=u, alpha=alpha, beta=beta)
end

# Finally, we define a function to plot the skyrmion velocity
function plot_velocity(us)
    vs = Float64[]
    js = Float64[]
    for u in us
        #for each u, we compute the velocity
        filename = @sprintf("assets/pos_u_%g_alpha_%g_beta_%g.txt", u, alpha, beta)
        data = readdlm(filename)
        ts, X, Y = data[:,1], data[:,2],  data[:,3]
        v = (last(X) - first(X))/(last(ts)-first(ts))
        push!(vs, v)
        push!(js, u*1.282)
    end
    print(vs, js)
    
    fig = Figure(resolution = (800, 480), fontsize=28)
    ax = Axis(fig[1, 1],
        xlabel = "j (10^10 A/m^2)",
        ylabel = "Velocity (m/s)"
    )

    #scatterlines!(ax, js, vs, label="X")
    scatterlines!(ax, js, vs, markersize = 8)
    #axislegend()
    save("assets/velocity.png",fig)
    return fig

end

plot_velocity(us)


