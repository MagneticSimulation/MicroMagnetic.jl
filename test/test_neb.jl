using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function test_convert()
    r, theta, psi = 2, pi / 3, pi / 2
    m0 = [r, theta, psi]
    cartesian = MicroMagnetic.spherical2cartesian(m0, 1)
    @test isapprox(cartesian[1], r * sin(theta) * cos(psi))
    @test isapprox(cartesian[2], r * sin(theta) * sin(psi))
    @test isapprox(cartesian[3], r * cos(theta))

    spherical = MicroMagnetic.cartesian2spherical(cartesian, 1)
    @test isapprox(spherical[1], r)
    @test isapprox(spherical[2], theta)
    @test isapprox(spherical[3], psi)
end

function test_inner_product()
    x1 = [1, 2.0, 3.0]
    x2 = [4.0, 5.0, 6.0]
    z = MicroMagnetic.inner_product(x1, x2, 1)
    @test isapprox(z[1], 4 + 10 + 18)
end

function test_slerp()
    t1 = pi / 2
    x1 = MicroMagnetic.slerp([1, 0, 0.0], [cos(t1), sin(t1), 0.0], 1 / 3, 1)
    @test isapprox(x1, [cos(t1 / 3), sin(t1 / 3), 0.0], atol=1e-7)

    t2 = 5 / 4 * pi
    t2_c = t2 - 2 * pi
    x2 = MicroMagnetic.slerp([1, 0, 0.0], [cos(t2), sin(t2), 0.0], 2 / 3, 1)
    @test isapprox(x2, [cos(t2_c * 2 / 3), sin(t2_c * 2 / 3), 0.0])
end

function test_init_images()
    t1 = pi / 2
    x1 = [1, 0, 0.0]
    x2 = [cos(t1), sin(t1), 0.0]
    x3 = [-1, 0.0, 0.0]
    images = MicroMagnetic.init_images((x1, x2, x3), (2, 1), 1)
    @test isapprox(images[:, 1], x1)
    @test isapprox(images[:, 2], [cos(t1 / 3), sin(t1 / 3), 0.0], atol=1e-7)
    @test isapprox(images[:, 3], [cos(2 * t1 / 3), sin(2 * t1 / 3), 0.0], atol=1e-7)
    @test isapprox(images[:, 4], [cos(t1), sin(t1), 0.0])
    @test isapprox(images[:, 5], [cos(3 * pi / 4), sin(3 * pi / 4), 0.0], atol=1e-7)
    @test isapprox(images[:, 6], [-1, 0.0, 0.0])
end

function test_init_neb()
    function create_sim()
        mesh = FDMesh(; nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
        sim = Sim(mesh; name="neb", driver="None", save_data=false)
        set_Ms(sim, 8e5)

        init_m0(sim, (0.6, 0, -0.8))

        add_exch(sim, 1.3e-11)
        add_anis(sim, 5e4; axis=(1, 0, 0), name="Kx")
        add_anis(sim, -2e4; axis=(0, 0, 1), name="Kz")
        return sim
    end
    sim = create_sim()
    init_m0(sim, (1, 0, 0.0))
    m0 = Array(sim.spin)

    init_m0(sim, (0, 1, 0.0))
    m1 = Array(sim.spin)
    neb = MicroMagnetic.NEB(sim, (m0, m1), (1,); name="test")

    @test isapprox(neb.spin, [0.7071067811865475, 0.7071067811865475, 0.0])
end

function test_neb_macrospin()
    function creat_sim()
        mesh = FDMesh(; nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
        sim = Sim(mesh; name="neb", driver="None", save_data=false)
        set_Ms(sim, 8e5)

        init_m0(sim, (0.6, 0, -0.8))

        #add_exch(sim, 1.3e-11)
        add_anis(sim, 5e4; axis=(1, 0, 0), name="Kx")
        add_anis(sim, -2e4; axis=(0, 0, 1), name="Kz")
        return sim
    end

    function init_m(i, j, k, dx, dy, dz)
        return (1, 0, 0)
    end

    sim = creat_sim()

    neb = NEB(sim, [(1, 0, 0), (0, 1, 1), (-1, 0, 0)], [5, 5]; name="test", driver="SD")
    @test(length(neb.spin) == 3 * (5 + 5 + 1))

    init_energy_diff = maximum(neb.energy_cpu) - neb.energy_cpu[1]
    relax(neb; max_steps=500, stopping_dmdt=0.001)
    final_energy_diff = maximum(neb.energy_cpu) - neb.energy_cpu[1]
    expected_energy_diff = 5e4 * 5e-9^3
    println("init_diff=$init_energy_diff, expected_diff=$expected_energy_diff, final_diff=$final_energy_diff")
    @test init_energy_diff > expected_energy_diff
    @test abs(final_energy_diff - expected_energy_diff) / expected_energy_diff < 1e-3
end

function test_neb()
    test_convert()
    test_inner_product()
    test_slerp()
    test_init_images()
    test_init_neb()
    return
end

@using_gpu()
test_functions("NEB", test_neb; platforms=["CPU"])
test_functions("NEB MacroSpin", test_neb_macrospin)
