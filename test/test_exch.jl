using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function A_fun(i, j, k, dx, dy, dz)
    return 1.3e-11
end

function test_exch_scalar(Nx=50)
    function m0_fun(i, j, k, dx, dy, dz)
        L = Nx * dx
        x = i * dx
        return sin(2 * pi * x / L), sin(2 * pi * x / L + 1.2), sin(2 * pi * x / L + 2.3)
    end

    mesh = FDMesh(; dx=2e-9, nx=Nx, ny=1, nz=1, pbc="x")

    Ms = 8.6e5
    A = 1.3e-11

    sim = Sim(mesh)
    set_Ms(sim, Ms)
    init_m0(sim, m0_fun; norm=false)
    exch = add_exch(sim, A_fun)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    if isa(sim.spin, Array{Float64})
        f1 = Array(exch.field)
        MicroMagnetic.effective_field_debug(exch, sim, sim.spin, 0.0)
        @test isapprox(f1, exch.field, atol=1e-10)
    end

    xs = (1:Nx) * 2e-9
    mu0 = 4 * pi * 1e-7
    L = Nx * 2e-9
    expected_x = -2 * A / (mu0 * Ms) * (2 * pi / L)^2 * sin.((2 * pi / L) .* xs)
    expected_y = -2 * A / (mu0 * Ms) * (2 * pi / L)^2 * sin.((2 * pi / L) .* xs .+ 1.2)
    expected_z = -2 * A / (mu0 * Ms) * (2 * pi / L)^2 * sin.((2 * pi / L) .* xs .+ 2.3)

    b = reshape(Array(exch.field), 3, sim.n_total)
    @test (maximum(b[1, :] .- expected_x) ./ Ms < 2e-4)
    @test (maximum(b[2, :] .- expected_y) ./ Ms < 2e-4)
    @test (maximum(b[3, :] .- expected_z) ./ Ms < 2e-4)

    Delta = 2e-9
    mesh = FDMesh(; dx=2e-9, dz=Delta, nx=1, ny=1, nz=3, pbc="x")
    sim = Sim(mesh)
    set_Ms(sim, Ms)

    sigma = 1e-5
    return init_m0(sim, (0.6, 0.8, 0))
end

function Ms_x(i, j, k, dx, dy, dz)
    return j == 2 && k == 2 ? 1e5 : 0
end

function Ms_y(i, j, k, dx, dy, dz)
    return i == 2 && k == 2 ? 1e5 : 0
end

function Ms_z(i, j, k, dx, dy, dz)
    return i == 2 && j == 2 ? 1e5 : 0
end

function test_exch_vector(direction=:x)
    mesh = FDMesh(; nx=3, ny=3, nz=3)
    function m_fun(i, j, k, dx, dy, dz)
        return (i^2, j + 1.0, k * j)
    end

    sim = Sim(mesh)

    if direction == :x
        set_Ms(sim, Ms_x)
        A = (1.3e-11, 0, 0)
    elseif direction == :y
        set_Ms(sim, Ms_y)
        A = (0, 1.3e-11, 0)
    elseif direction == :z
        set_Ms(sim, Ms_z)
        A = (0, 0, 1.3e-11)
    end

    init_m0(sim, m_fun)

    ex1 = add_exch(sim, A_fun; name="ex1")
    ex2 = add_exch(sim, A; name="ex2")

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    f1 = Array(ex1.field)
    f2 = Array(ex2.field)

    #print(maximum(abs.(f1.-f2)))
    @test isapprox(f1, f2, atol=1e-7)
end

function test_exch_vectors()
    test_exch_vector(:x)
    test_exch_vector(:y)
    return test_exch_vector(:z)
end

@using_gpu()
test_functions("Exchange", test_exch_scalar, test_exch_vectors)
