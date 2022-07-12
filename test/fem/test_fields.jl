using JuMag
using Test

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

test_zeeman()