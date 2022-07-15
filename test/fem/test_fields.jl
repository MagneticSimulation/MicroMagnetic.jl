using JuMag
using Test
using DelimitedFiles

function test_zeeman()
    mesh = FEMesh("octa.neutral")
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
    mesh = FEMesh("fields/cylinder.neutral")
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)
    N = length(sim.spin)
    println("The length of spin: ", N)
    
    m0 = readdlm("fields/m.txt", header=true)[1]
    mxyz = reshape(transpose(m0[:, 4:6]), N, 1)
    #println("max diff: ", maximum(abs.(sim.spin - mxyz)))
    eps = 1e-6
    @test maximum(abs.(sim.spin - mxyz)) < eps
    #print(mxyz[1:50])
end

function test_exchange_field()
    mesh = FEMesh("fields/cylinder.neutral")
    sim = Sim(mesh)
    set_Ms(sim, 8.6e5)

    init_m0(sim, init_m_fun)

    z = add_exch(sim, 1.3e-11)

    JuMag.effective_field(sim, sim.spin, 0.0)
    N = 3*sim.n_nodes

    f0 = readdlm("fields/exch.txt", header=true)[1]
    field = reshape(transpose(f0[:, 4:6]), N, 1)
    eps = 5e-5
    @test maximum(abs.(z.field - field)) < eps

end


test_zeeman()
test_init_m_function()
test_exchange_field()