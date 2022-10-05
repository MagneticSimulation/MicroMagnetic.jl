using JuMag
using Printf
using NPZ
using Test
JuMag.cuda_using_double(true)
function test_neb_gpu()
    m1=[1.0,0,0]
    m2=[0,1.0,0]
    m3=[0,0,1.0]

    mxy=JuMag.rotate_to(m1,m2,pi/4)
    mxz=JuMag.rotate_to(m1,m3,pi/4)
    myz=JuMag.rotate_to(m2,m3,pi/4)

    function test_rotation_operator()
        @test isapprox(JuMag.rotation_operator(m1,m3,pi/2),m2)
        @test isapprox(JuMag.rotation_operator(m1,m3,pi/4),[sqrt(2)/2,sqrt(2)/2,0])
        @test isapprox(JuMag.rotation_operator(m3,m1,pi/2),-m2)
    end
    function test_rotate_to()
        @test isapprox(mxy,[sqrt(2)/2,sqrt(2)/2,0])
        @test isapprox(mxz,[sqrt(2)/2,0,sqrt(2)/2])
        @test isapprox(myz,[0,sqrt(2)/2,sqrt(2)/2])
    end
    test_rotation_operator()
    @info("test_rotation_operator passed!")
    test_rotate_to()
    @info("test_rotate_to passed!")
    function test_interpolation()
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
        m_neb = JuMag.interpolate_m(m1, m2, 1)
        println("Input: ", m1)
        println("interpolate_m1: ", m_neb)
        p = JuMag.interpolate_m_spherical(m1, m2, 1)
        println("interpolate_m2: ", p)
    end
    function init_m(i,j,k,dx,dy,dz)
        return (1,0,0)
    end
    function creat_sim()
        mesh =  FDMeshGPU(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
        sim = Sim(mesh, name="neb", driver="None", save_data=false)
        set_Ms(sim, 8e5)

        init_m0(sim, (0.6, 0, -0.8))

        add_exch(sim, 1.3e-11)
        add_anis(sim, 5e4, axis=(1,0,0), name="Kx")
        add_anis(sim, -2e4, axis=(0,0,1), name="Kz")
        return sim
    end

    sim=creat_sim()

    neb = NEB_GPU(sim, [ (1,0,0), (0,1,0)], 1; name="test")
    println("left:",neb.image_l, typeof(neb.image_l))
    println("right:",neb.image_r)
    println(neb.spin)
    println("distance:",neb.distance)
    @test isapprox(neb.distance[1],atan(sqrt(2)/2,sqrt(2)/2))
    @test isapprox(Array(neb.image_l), [1,0,0])
    @test isapprox(Array(neb.spin[1:3]), mxy)
    @test isapprox(Array(neb.image_r),[0,1,0])

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
