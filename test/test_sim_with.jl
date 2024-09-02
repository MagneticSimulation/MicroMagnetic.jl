using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function analytical_llg(alpha::Float64, gamma::Float64, H0::Float64, t::Float64)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * t

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return [mx, my, mz]
end

function test_sim_with()
    mesh = FDMesh(; nx=2, ny=2, nz=2)
    args = (task="Relax", mesh=mesh, m0=(1, 2, 3), Ms=8e5, A=1e-11, H=(0, 0, 1e5),
            name="test_relax", stopping_dmdt=0.01)

    sim = sim_with(args)
    m = Array(sim.spin)
    @test isapprox(m[1:3], [0, 0, 1])

    mesh = FDMesh(; nx=1, ny=1, nz=1)
    analytical = analytical_llg(0.05, 2.21e5, 1e5, 3e-10)
    args = (task="Dynamics", mesh=mesh, m0=(1, 0, 0), Ms=8e5, H=(0, 0, 1e5), gamma=2.21e5,
            alpha=0.05, steps=10, dt=3e-11, tol=1e-8)

    sim = sim_with(args)
    m = Array(sim.spin)
    @test isapprox(m[1:3], analytical)

    args = (task_s=["relax", "dynamics"], mesh=mesh, m0=(1, 2, 3), Ms=8e5, A=1e-11, steps=10,
            alpha=0.05, tol=1e-8, dt=3e-11, stopping_dmdt=0.001,
            H_s=[(1e5, 0, 0), (0, 0, 1e5)])

    sim = sim_with(args)
    m = Array(sim.spin)
    println(m[1:3])
    @test isapprox(m[1:3], analytical)

    args = (name="test_read", task="relax", mesh=mesh, m0=(1, 2, 3), Ms=8e5, A=1e-11, stopping_dmdt=0.01,
            H_s=[(1e5, 0, 0), (0, 0, -1e5)])

    sim = sim_with(args)
    m = Array(sim.spin)
    @test isapprox(m[1:3], [0, 0, -1])
    return nothing
end

function test_read_table()
        data, units = read_table("test_read_sd.txt")
        @test data["m_x"][end] < 1e-15
        @test data["m_y"][end] < 1e-15
        @test data["m_z"][end]+1 < 1e-15
end


#@using_gpu()
test_functions("sim_with", test_sim_with)
test_functions("read_table", test_read_table)