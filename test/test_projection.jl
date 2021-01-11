using JuMag
using Test
using Printf


function test_tilt_beta()
	mesh = FDMesh(nx=128,ny=128,nz=1)
	sim = Sim(mesh)
	set_Ms(sim,1e4)
	init_m0(sim,(0,0,1))
	m = sim.spin
	smx,smy,smz = JuMag.Make_Projection(m,128,128,1,Nx=128,Ny=128,beta=pi/3,gamma=0)
	@test(abs(smx[64,64]+sqrt(3)/2*2)<1e-4)
	@test(abs(smx[56,56]+sqrt(3)/2*2)<1e-4)
	@test(abs(smz[64,64]-1/2*2)<1e-4)
	@test(abs(smz[56,56]-1/2*2)<1e-4)
	@test(abs(sum(smx)+sqrt(3)/2*128*128)<1e-2)
end

function test_tilt_gamma()
	mesh = FDMesh(nx=128,ny=128,nz=1)
	sim = Sim(mesh)
	set_Ms(sim,1e4)
	init_m0(sim,(0,0,1))
	m = sim.spin
	smx,smy,smz = JuMag.Make_Projection(m,128,128,1,Nx=128,Ny=128,beta=0,gamma=pi/3)
	@test(abs(smy[64,64]+sqrt(3)/2*2)<1e-4)
	@test(abs(smy[56,56]+sqrt(3)/2*2)<1e-4)
	@test(abs(smz[64,64]-1/2*2)<1e-4)
	@test(abs(smz[56,56]-1/2*2)<1e-4)
	@test(abs(sum(smy)+sqrt(3)/2*128*128)<1e-2)
end
function test_tiltbeta_multilayer()
	mesh = FDMesh(nx=64,ny=64,nz=64)
	sim = Sim(mesh)
	set_Ms(sim,1e4)
	init_m0(sim,(0,0,1))
	m = sim.spin
	smx,smy,smz = Make_Projection(m,64,64,64,Nx=128,Ny=128,beta=pi/4,gamma=0)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64+smx[64,64]) )<1e-1)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64+smx[64,50]) )<1e-1)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64-smz[64,64]) )<1e-1)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64-smz[64,50]) )<1e-1)
	@test(sum(smx)+64*64*64*sqrt(2)/2<1)
end
function test_tiltgamma_multilayer()
	mesh = FDMesh(nx=64,ny=64,nz=64)
	sim = Sim(mesh)
	set_Ms(sim,1e4)
	init_m0(sim,(0,0,1))
	m = sim.spin
	smx,smy,smz = Make_Projection(m,64,64,64,Nx=128,Ny=128,gamma=pi/4)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64+smy[64,64]) )<1e-1)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64+smy[50,64]) )<1e-1)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64-smz[64,64]) )<1e-1)
	@test(abs( (32*sqrt(2)-0.5)/(32*sqrt(2)*64-smz[50,64]) )<1e-1)
	@test(sum(smy)+64*64*64*sqrt(2)/2<1)
end

@testset "tools" begin
	test_Normal_Projection()
	mini_test2()
	test_tilt_beta()
	test_tilt_gamma()
	test_tiltbeta_multilayer()
	test_tiltgamma_multilayer()
end