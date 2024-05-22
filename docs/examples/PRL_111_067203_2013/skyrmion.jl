# ---
# title: Magnetic skyrmion
# author: Weiwei Wang
# date: 2022-10-04
# description: to check how to use the parameters given in PRL **111** 067203 (2013).
# tag: atomistic; skyrmion
# ---

#==============
!!! note "Used parameters in the simulation"
    |Parameter | Value  | 
    | :----:   | :----: | 
    | Lattice constant | $a = 0.5$ nm |
    | Spin length      | $S = 1$      | 
    | Magnetic moment  |  $\mu_s = 2 \mu_B$ |
    | Excahnge constant |  $J = 50 k_B$   |
    | DMI         | $D/J = 0.5$  |     |
    | External field  | $H \mu_s /J  = 0.2$ | 

===============#

using MicroMagnetic
using Printf
using NPZ


mesh =  CubicMesh(nx=150, ny=50, nz=1, pbc="xy")

function m0_fun(i,j,k, dx, dy, dz)
  r2 = (i-70)^2 + (j-25)^2
  if r2 < 10^2
    return (0.01, 0, -1)
  end
  return (0,0,1)
end

function relax_system()
  sim = Sim(mesh, driver="SD", name="skx")
  set_mu_s(sim, mu_s_1)
  init_m0(sim, m0_fun)

  J = 50 * k_B
  add_exch(sim, J, name="exch")
  add_dmi(sim, 0.5*J, name="dmi")
  
  Hz = 0.2 * J / mu_s_1
  add_zeeman(sim, (0,0,Hz))
  
  relax(sim, maxsteps=2000, stopping_dmdt=0.01)

  save_vtk(sim, "skx", fields=["exch", "dmi"])
  
end

relax_system()
