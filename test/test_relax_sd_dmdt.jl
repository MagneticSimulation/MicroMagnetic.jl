using MicroMagnetic
using Test
using Statistics

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

# The SD driver exposes the raw spherical-gradient norm max|gk| of the state it
# just stepped from (driver.max_gk, pinned spins masked); relax scales it with
# the gyromagnetic ratio of the sim's field convention (_sd_gamma: A/m fields ->
# gamma0/2, Tesla fields -> gamma0/(2*mu_0)) into an LLG-comparable dmdt. These
# tests pin the calibration for both conventions on a single-spin problem where
# |gk| = |H_zeeman| exactly (m starts perpendicular to the field, and the
# anisotropy/exchange fields vanish or stay parallel to m).

# one SD step's relax-style dmdt reading, in 1/s of the sim's convention
function sd_first_step_dmdt(sim)
    MicroMagnetic.run_step(sim, sim.driver)
    return sim.driver.max_gk * MicroMagnetic._sd_gamma(sim) / 2
end

function test_sd_dmdt_micromagnetic()
    # m = (1, 0.2, 0) has mz = 0, so anisotropy contributes nothing and |gk| = H_z
    mesh = FDMesh(; nx=4, ny=4, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh; driver="SD", name="sd_dmdt_micro", save_data=false, quiet=true)
    set_Ms(sim, 4.2e5)
    add_exch(sim, 10e-12)
    add_anis(sim, 1.61e5; axis=(0, 0, 1))
    add_zeeman(sim, (0, 0, -200 * mT))
    init_m0(sim, (1, 0.2, 0))

    dmdt = sd_first_step_dmdt(sim)
    expected = (0.2 / mu_0) * 2.21e5 / 2      # |H| in A/m times gamma0/2
    @test isapprox(dmdt, expected; rtol=0.05)
end

function test_sd_dmdt_atomistic()
    # same geometry in the atomistic convention: field in Tesla, gamma0/mu_0
    mesh = CubicMesh(; nx=2, ny=2, nz=1, dx=1e-9)
    sim = Sim(mesh; driver="SD", name="sd_dmdt_at", save_data=false, quiet=true)
    set_mu_s(sim, 2 * mu_B)
    add_exch(sim, 1.0e-20)
    add_anis(sim, 1.0e-22; axis=(0, 0, 1))
    add_zeeman(sim, (0, 0, -0.2))             # Tesla
    init_m0(sim, (1, 0.2, 0))

    dmdt = sd_first_step_dmdt(sim)
    expected = 0.2 * (2.21e5 / mu_0) / 2      # |B| in T times gamma0/mu_0/2
    @test isapprox(dmdt, expected; rtol=0.05)

    # end to end: a 0.05 deg/ns threshold must not trigger during the large
    # rotations, and the criterion run must land on the unconstrained
    # relaxation's state
    function mk(name)
        s = Sim(mesh; driver="SD", name=name, save_data=false, quiet=true)
        set_mu_s(s, 2 * mu_B)
        add_exch(s, 1.0e-20)
        add_anis(s, 1.0e-22; axis=(0, 0, 1))
        add_zeeman(s, (0, 0, -0.2))
        init_m0(s, (1, 0.2, 0))
        return s
    end
    crit = mk("sd_dmdt_crit")
    relax(crit; max_steps=2000, stopping_dmdt=0.05)
    free = mk("sd_dmdt_free")
    relax(free; max_steps=2000, stopping_dmdt=-1.0)

    @test crit.driver.steps > 10
    # the anisotropy wells at mz = ±1 are nearly degenerate (the 0.2 T field only
    # tilts them), so which one SD lands in depends on rounding noise — assert a
    # converged single-domain state and agreement with the unconstrained run
    @test maximum(abs.(Array(crit.spin) - Array(free.spin))) < 1e-5
    @test abs(abs(mean(Float64.(Array(crit.spin))[3:3:end])) - 1) < 1e-5
end

function test_relax_sd_dmdt()
    test_sd_dmdt_micromagnetic()
    test_sd_dmdt_atomistic()
    return nothing
end

@using_gpu()
test_functions("RelaxSDdmdt", test_relax_sd_dmdt)
