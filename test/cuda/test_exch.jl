using JuMag
using Test

Nx = 50

function m0_fun(i,j,k,dx,dy,dz)
  L = Nx*dx
  x = i*dx
  return sin(2*pi*x/L), sin(2*pi*x/L+1.2), sin(2*pi*x/L+2.3)
end

#Test mesh
mesh =  FDMeshGPU(dx=2e-9, nx=Nx, ny=1, nz=1, pbc="x")
@test isapprox(mesh.dx, 2e-9, rtol=1e-7)
@test mesh.nx == Nx

Ms = 8.6e5
A = 1.3e-11

sim = Sim(mesh)
set_Ms(sim, Ms)

init_m0(sim, m0_fun, norm=false)
add_exch(sim, A)

JuMag.effective_field(sim, sim.spin, 0.0)

xs = (1:Nx)*2e-9
mu0 = 4*pi*1e-7
L = Nx*2e-9
expected_x =-2*A/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs);
expected_y =-2*A/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs.+1.2);
expected_z =-2*A/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs.+2.3);

field = Array(sim.field)
println(field)
println(expected_x)
b = reshape(field, 3, sim.nxyz)
println(maximum(b[1,:].-expected_x)./Ms)
@test (maximum(b[1,:].-expected_x)./Ms<2e-4)
@test (maximum(b[2,:].-expected_y)./Ms<2e-4)
@test (maximum(b[3,:].-expected_z)./Ms<2e-4)

#test exchange anistropy
kea = 1e-11
sim = Sim(mesh)
set_Ms(sim, Ms)

init_m0(sim, m0_fun, norm=false)
add_exch_anis(sim, kea)

JuMag.effective_field(sim, sim.spin, 0.0)

xs = (1:Nx)*2e-9
mu0 = 4*pi*1e-7
L = Nx*2e-9

expected_x =-2*kea/(mu0*Ms)*(2*pi/L)^2*sin.((2*pi/L).*xs);
expected_y =0.0;
expected_z =0.0;

field = Array(sim.field)
println(field)
println(expected_x)
b = reshape(field, 3, sim.nxyz)
println(maximum(b[1,:].-expected_x)./Ms)
@test (maximum(b[1,:].-expected_x)./Ms<2e-4)
@test (maximum(b[2,:].-expected_y)./Ms<2e-4)
@test (maximum(b[3,:].-expected_z)./Ms<2e-4)

# test rkky
JuMag.cuda_using_double(true)
Delta = 2e-9
mesh =  FDMeshGPU(dx=2e-9, dz=Delta, nx=1, ny=1, nz=3, pbc="x")
sim = Sim(mesh)
set_Ms(sim, Ms)

sigma = 1e-5

function m0_fun_rkky(i,j,k,dx,dy,dz)
  if k == 1
    return 0.6, 0.8, 0
  elseif k == 3
    return -0.8, 0, 0.6
  else
    return 1,0,0
  end
end

init_m0(sim, m0_fun_rkky)
r = add_exch_rkky(sim, sigma)

#JuMag.effective_field(sim, sim.spin, 0.0)
JuMag.compute_system_energy(sim, sim.spin, 0.0)
b = Array(sim.field)
println("field = ", b)
f1 = sigma/Delta/(mu0*Ms)*0.6
f2 = sigma/Delta/(mu0*Ms)*0.8
expected = [-f2, 0, f1, 0,0,0, f1, f2, 0]
println(expected-b)
println("max diff: ", maximum(abs.(expected-b)))
@test maximum(abs.(expected-b)) < eps()
println("energy: ",r.total_energy)
