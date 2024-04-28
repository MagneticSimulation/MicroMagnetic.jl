# ---
# title: Magnetic skyrmion
# author: Weiwei Wang
# date: 2022-10-05
# description: to check how to use the parameters given in Nat Commun **4** 1463 (2013).
# tag: atomistic; skyrmion
# ---

#==============
In the work, the authors considered the classical Heisenberg model on the 2D square lattice.
The local magnetic moment $\mathbf{M}$ is the local spin $\mathbf{S}$ in JuMag. 

!!! note "Used parameters in the simulation"
    |Parameter | Value  | 
    | :----:   | :----: | 
    | Lattice constant | $a = 0.5$ nm |
    | Spin length      | $S = 1$      | 
    | Magnetic moment  |  $\mu_s = 2 \mu_B$ |
    | Excahnge constant |  $J = 1$ meV   |
    | DMI         | $D/J = 0.18$  |     |
    | External field  | $H \mu_s /J  = 0.015$ | 
    | Spin-polarization rate        |  $p = 0.2$   | 

===============#

using JuMag
using Printf
using NPZ

function m0_fun(i, j, k, dx, dy, dz)
    r = 25
    if ((i - 80)^2 + (j - 40)^2 < r^2)
        return (0.05, 0.01, -1)
    end
    return (0,0,1)
end

function relax_system()
    J = 1*meV
    D = 0.18*J

    mesh = CubicMesh(nx=166, ny=80, nz=1, dx=0.5e-9, dy=0.5e-9, pbc="xy")

    sim = Sim(mesh, driver="SD", name="skx")
    set_mu_s(sim, mu_s_1)
    
    add_exch(sim, J, name="exch")
    add_dmi(sim, D, name="dmi")
    
    Hz= 1.5e-2*J / mu_s_1
    add_zeeman(sim, (0,0,Hz))
    init_m0(sim, m0_fun)

    relax(sim, maxsteps=1000, stopping_dmdt=1e-5, using_time_factor=false)

    save_vtk(sim, "skx") 
    
    Q = compute_skyrmion_number(Array(sim.spin),mesh)
  
    return Q
end

relax_system()
