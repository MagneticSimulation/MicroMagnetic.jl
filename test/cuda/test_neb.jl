using JuMag
using Printf
using NPZ
using Test
JuMag.cuda_using_double(true)
function test_neb_gpu()
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
    function init_m(i,j,k,dx,dy,dz)
        return (1,0,0)
    end
    function creat_sim()
        mesh =  FDMeshGPU(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
        sim = Sim(mesh, name="neb", driver="none",save_data=false)
        set_Ms(sim, 8e5)

        init_m0(sim, (0.6, 0, -0.8))

        add_exch(sim, 1.3e-11)
        add_anis(sim, 5e4, axis=(1,0,0), name="Kx")
        add_anis(sim, -2e4, axis=(0,0,1), name="Kz")
        return sim
    end

    sim=creat_sim()

    neb = NEB_GPU(sim, [ (1,0,0), (0,1,0)], 1; name="test")
    println("left:",neb.image_l)
    println("right:",neb.image_r)
    println(neb.spin)
    println("distance:",neb.distance)
    @test isapprox(neb.distance[1],atan(sqrt(2)/2,sqrt(2)/2))
    @test isapprox(neb.image_l[:],[1,0,0])
    @test isapprox(neb.spin[1:3],mxy)
    @test isapprox(neb.image_r[:],[0,1,0])

    neb = NEB_GPU(sim, [ (1,0,0), (0,1,1), (-1,0,0)], [5, 5]; name="test")
    @test(size(neb.spin) == (33,))

    init_energy_diff = maximum(neb.energy_cpu) - neb.energy_cpu[1]
    relax(neb, maxsteps=200, stopping_dmdt=0.01, save_ovf_every=-1)
    final_energy_diff = maximum(neb.energy_cpu) - neb.energy_cpu[1]
    expected_energy_diff = 5e4*5e-9^3
    println("init_diff=$init_energy_diff, expected_diff=$expected_energy_diff, final_diff=$final_energy_diff")
    @test init_energy_diff > expected_energy_diff
    @test abs(final_energy_diff - expected_energy_diff)/expected_energy_diff < 1e-8
end
test_neb_gpu()