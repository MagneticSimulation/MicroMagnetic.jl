using JuMag
using Printf
using NPZ
using Test

m1=[1.0,0,0]
m2=[0,1.0,0]
m3=[0,0,1.0]

mxy=JuMag.rotation_operator(m1,m2,pi/4)
mxz=JuMag.rotation_operator(m1,m3,pi/4)
myz=JuMag.rotation_operator(m2,m3,pi/4)
@test isapprox(mxy,[sqrt(2)/2,sqrt(2)/2,0])
@test isapprox(mxz,[sqrt(2)/2,0,sqrt(2)/2])
@test isapprox(myz,[0,sqrt(2)/2,sqrt(2)/2])
interplotion1=JuMag.interpolate_m([1.0,0,0],[0,1.0,0],1)
interplotion2=JuMag.interpolate_m([1.0,0,0],[1.0,0,0],1)
interplotion3=JuMag.interpolate_m([0,0,1.0],[0,0,1.0],1)
@test isapprox(interplotion1[:,1],[1.0,0,0])
@test isapprox(interplotion1[:,2],[sqrt(2)/2,sqrt(2)/2,0])
@test isapprox(interplotion2[:,1],[1.0,0,0])
@test isapprox(interplotion2[:,2],[1.0,0,0])
@test isapprox(interplotion3[:,1],[0,0,1.0])
@test isapprox(interplotion3[:,2],[0,0,1.0])

m1 = [0.061446725558210215, 0.32475341620029013, 0.9438005714050057]
m2 = [0, 0, 1.0]
m = JuMag.interpolate_m(m1, m2, 1)
println("Input: ", m1)
println("interpolate_m1: ", m)
p = JuMag.interpolate_m_spherical(m1, m2, 1)
println("interpolate_m2: ", p)

function creat_sim()
    mesh =  FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh, name="neb", driver="none",save_data=false)
    set_Ms(sim, 8e5)

    init_m0(sim, (0.6, 0, -0.8))

    add_exch(sim, 1.3e-11)
    add_anis(sim, 5e4, axis=(1,0,0), name="Kx")
    add_anis(sim, -2e4, axis=(0,0,1), name="Kz")
    return sim
end

function init_m(i,j,k,dx,dy,dz)
    return (1,0,0)
end

sim=creat_sim()

neb = NEB(sim, [ (1,0,0), (0,1,0)], 1; name="test")
println(neb.images)
@test isapprox(neb.distance[1],atan(sqrt(2)/2,sqrt(2)/2))
@test isapprox(neb.images[:,1],[1,0,0])
@test isapprox(neb.images[:,2],mxy)
@test isapprox(neb.images[:,3],[0,1,0])

neb = NEB(sim, [ (1,0,0), (0,1,1), (-1,0,0)], [5, 5]; name="test")
@test(size(neb.images) == (3, 13))

init_energy_diff = maximum(neb.energy) - neb.energy[1]
relax(neb, maxsteps=200, stopping_dmdt=0.01, save_ovf_every=-1)
final_energy_diff = maximum(neb.energy) - neb.energy[1]
expected_energy_diff = 5e4*5e-9^3
println("init_diff=$init_energy_diff, expected_diff=$expected_energy_diff, final_diff=$final_energy_diff")
@test init_energy_diff > expected_energy_diff
@test abs(final_energy_diff - expected_energy_diff)/expected_energy_diff < 1e-8
