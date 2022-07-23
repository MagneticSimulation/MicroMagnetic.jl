using JuMag
using Test
using DelimitedFiles

function test_zeeman()
    filepath = joinpath(@__DIR__, "octa.neutral")
    mesh = FEMesh(filepath)
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)
    
    z = add_zeeman(sim, (1,2,1e5))

    JuMag.effective_field(sim, sim.spin, 1.23e-11)

    @test z.field[1] == 1.0
    @test z.field[2] == 2.0
    @test z.field[3] == 1e5
end

function init_m_fun(x,y,z)
    r = sqrt(x*x+y*y+z*z)
    return (0, sin(r), cos(r))
end

function test_init_m_function()
    filepath = joinpath(@__DIR__, "fields/cylinder.neutral")
    mesh = FEMesh(filepath)
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)
    N = length(sim.spin)
    println("The length of spin: ", N)
    
    filepath = joinpath(@__DIR__, "fields/m.txt")
    m0 = readdlm(filepath, header=true)[1]
    mxyz = reshape(transpose(m0[:, 4:6]), N, 1)
    #println("max diff: ", maximum(abs.(sim.spin - mxyz)))
    eps = 1e-6
    @test maximum(abs.(sim.spin - mxyz)) < eps
    #print(mxyz[1:50])
end


function test_zeeman_field()
    filepath = joinpath(@__DIR__, "fields/cylinder.neutral")
    mesh = FEMesh(filepath)
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_zeeman(sim, (0,0,1e5))

    JuMag.effective_field(sim, sim.spin, 0.0)
    
    expected_energy = 1.81151949357e-21
    println("zeeman energy: ",sum(z.energy))
    @test abs(sum(z.energy)-expected_energy)/expected_energy<1e-10
end

function test_anis_field()
    filepath = joinpath(@__DIR__, "fields/cylinder.neutral")
    mesh = FEMesh(filepath)
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_anis(sim, 5.2e5, axis=(0,0.6,0.8))

    JuMag.effective_field(sim, sim.spin, 0.0)
    N = 3*sim.n_nodes

    filepath = joinpath(@__DIR__, "fields/anis.txt")
    f0 = readdlm(filepath, header=true)[1]
    field = reshape(transpose(f0[:, 4:6]), N, 1)
    eps = 1e-6
    @test maximum(abs.(z.field - field)) < eps

    #expected_energy = 9.60333028541e-21

    #println("anis energy: ",sum(z.energy))
    #@test (sum(z.energy)-expected_energy)/expected_energy<1e-10

end



function test_exchange_field()
    filepath = joinpath(@__DIR__, "fields/cylinder.neutral")
    mesh = FEMesh(filepath)
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_exch(sim, 1.3e-11)

    JuMag.effective_field(sim, sim.spin, 0.0)
    N = 3*sim.n_nodes

    filepath = joinpath(@__DIR__, "fields/exch.txt")
    f0 = readdlm(filepath, header=true)[1]
    field = reshape(transpose(f0[:, 4:6]), N, 1)
    eps = 5e-5
    @test maximum(abs.(z.field - field)) < eps

    expected_energy = 3.18671455452e-19

    println("exch energy: ",sum(z.energy))
    @test (sum(z.energy)-expected_energy)/expected_energy<1e-10

end

function test_demag_field()
    filepath = joinpath(@__DIR__, "fields/cylinder.neutral")
    mesh = FEMesh(filepath)
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    demag = add_demag(sim)

    JuMag.effective_field(sim, sim.spin, 0.0)
    N = 3*sim.n_nodes

    filepath = joinpath(@__DIR__, "fields/g_phi.txt")
    g_phi = readdlm(filepath, header=true)[1]
    g1 = g_phi[:, 4]
    phi1 = g_phi[:, 5]
    println("diff:", g1[1:10], demag.g1[1:10])
    println("g1 max diff: ", maximum(abs.(demag.g1 - g1)))
    @test maximum(abs.(demag.g1 - g1)) < 1e-9

    println("phi1 max diff: ", maximum(abs.(demag.phi1 - phi1))) 
    @test maximum(abs.(demag.phi1 - phi1)) < 5e-9

    filepath = joinpath(@__DIR__, "fields/demag.txt")
    f0 = readdlm(filepath, header=true)[1]
    field = reshape(transpose(f0[:, 4:6]), N, 1)
    eps = 5e-5
    mean = sum(abs.(field))/length(field)
    @test maximum(abs.(demag.field - field)) / mean < eps

    energy = 3.97e-21
    @test (sum(demag.energy)-energy)/energy<0.02

end

test_zeeman()
test_init_m_function()
test_zeeman_field()
test_anis_field()
test_exchange_field()
test_demag_field()