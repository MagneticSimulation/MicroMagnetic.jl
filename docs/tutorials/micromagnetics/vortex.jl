# ---
# title: Magnetic vortex
# author: Weiwei Wang
# date: 2022-05-19
# description: an example to demostrate how to obtain the vortex in JuMag.
# tag: tutorial
# ---

# We import JuMag and use double float precision in the simulation
using JuMag
using CairoMakie

# In this example, the studied system is a round nanodisk. Since we are using finite different method,
# we create a FDMesh with dimensions 200nm x 200nm x 20 nm. 
mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=5e-9, nx=100, ny=100, nz=4);

# Create a cylindrical shape with a diameter of 100 nm
geo = Cylinder(radius=100e-9); 

# To start the simulation, we need to give an initial state.
# Here we use a function to give the initial state. The function should take six parameters `(i,j,k,dx,dy,dz)`.
function init_fun(i,j,k,dx,dy,dz)
  x = i-50.5
  y = j-50.5
  r = (x^2+y^2)^0.5
  if r<20
    return (0,0,1)
  end
  return (y/r, -x/r, 0)
end

# Create a Sim instance and then relax the system
sim = create_sim(mesh, shape=geo, Ms=8e5, A = 1.3e-11, demag=true, m0=init_fun);
relax(sim, maxsteps=5000, stopping_dmdt=0.1)
save_vtk(sim, "vortex")

# We plot the magnetization distribution using the plot_m function
plot_m(sim, figsize=(400, 400), arrows=(30,30))