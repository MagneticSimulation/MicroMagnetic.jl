using Test
using MicroMagnetic

@testset "Global state guards" begin
    # pin the baseline so the assertions below don't depend on test order
    saved_n_sims = MicroMagnetic._n_sims[]
    MicroMagnetic._n_sims[] = 0
    set_precision(Float64)
    set_backend("cpu")

    mesh = FDMesh(nx=4, ny=4, nz=1)
    sim = Sim(mesh)
    @test eltype(sim.spin) == Float64

    # once a Sim exists, set_precision/set_backend warn that only new Sims are affected
    @test_warn "existing Sim" set_precision(Float32)
    @test_warn "existing Sim" set_backend("cpu")

    # the existing Sim keeps its precision; Sims created afterwards pick the new one up
    @test eltype(sim.spin) == Float64
    sim32 = Sim(mesh; name="gs_f32")
    @test eltype(sim32.spin) == Float32

    # only Float64/Float32/AbstractFloat are accepted
    @test_throws ErrorException set_precision(Float16)
    @test_throws ErrorException set_precision(BigFloat)

    set_precision(Float64)
    MicroMagnetic._n_sims[] = saved_n_sims
end
