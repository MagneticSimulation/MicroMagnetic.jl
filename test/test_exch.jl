using JuMag
using Test

Nx = 50

function m0_fun(i,j,k,dx,dy,dz)
  L = Nx*dx
  x = i*dx
  return sin(2*pi*x/L), sin(2*pi*x/L+1.2), sin(2*pi*x/L+2.3)
end

#Test mesh
mesh =  FDMesh(dx=2e-9, nx=Nx, ny=1, nz=1, pbc="x")
@test mesh.dx == 2e-9
@test mesh.nx == Nx
println(mesh.ngbs)
@test mesh.ngbs[1,1] == Nx
@test mesh.ngbs[1,Nx] == Nx-1
@test mesh.ngbs[2,1] == 2
@test mesh.ngbs[2,Nx] == 1

Ms = 8.6e5
A = 1.3e-11

sim = Sim(mesh)
sim.Ms[:] .= Ms

init_m0(sim, m0_fun, norm=false)
add_exch(sim, A)

JuMag.effective_field(sim, sim.spin, 0.0)

xs = (1:Nx)*2e-9
mu0 = 4*pi*1e-7
L = Nx*2e-9
expected_x =-2*A/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs);
expected_y =-2*A/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs.+1.2);
expected_z =-2*A/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs.+2.3);

b = reshape(sim.field, 3, sim.nxyz)
println(maximum(b[1,:].-expected_x)./Ms)
@test (maximum(b[1,:].-expected_x)./Ms<2e-4)
@test (maximum(b[2,:].-expected_y)./Ms<2e-4)
@test (maximum(b[3,:].-expected_z)./Ms<2e-4)
