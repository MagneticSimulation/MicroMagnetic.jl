using JuMag
using Test

#Test mesh
function test_anis(;gpu=false)
    if gpu
        mesh =  FDMeshGPU(nx=10, ny=1, nz=1)
    else
        mesh =  FDMesh(nx=10, ny=1, nz=1)
    end

    sim = Sim(mesh)
    Ms = 8.6e5
    Kc = 1e3
    set_Ms(sim, Ms)
    init_m0(sim, (0.6,0.8,0), norm=false)

    anis = add_cubic_anis(sim, Kc)
    if gpu
        JuMag.compute_fields_to_gpu(sim, sim.spin, 0.0)
    else
        JuMag.effective_field(sim, sim.spin, 0.0)
    end

    @test isapprox(anis.field[1], 1/(JuMag.mu_0*Ms)*4*Kc*0.6^3)
    @test isapprox(anis.field[2], 1/(JuMag.mu_0*Ms)*4*Kc*0.8^3)
    @test isapprox(anis.field[3], 1/(JuMag.mu_0*Ms)*4*Kc*0^3)
    @test isapprox(anis.energy[1], -Kc*(0.6^4+0.8^4)*1e-27)
end
