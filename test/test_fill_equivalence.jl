using MicroMagnetic
using Test

# The correctness main evidence of the uniform-parameter representation. For the
# same driver and grid, a uniform parameter stored as an O(1) `Fill` and the same
# values materialised into a dense array must produce *bit-identical* trajectories:
# both go through the same kernel, and the Fill's getindex folds into the same
# immediate as the dense load.

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function run_case(field::Symbol, uniform::Bool)
    mesh = FDMesh(nx=16, ny=16, nz=2, dx=3e-9, dy=3e-9, dz=3e-9)
    name = "equiv_$(field)_$(uniform ? "fill" : "dense")"
    sim = Sim(mesh; name=name, driver="LLG")
    n = sim.n_total

    if field == :alpha
        set_Ms(sim, 8e5)
        add_exch(sim, 1.3e-11)
        add_zeeman(sim, (0, 0, 1e5))
        uniform ? set_alpha(sim, 0.1) : set_alpha(sim, fill(0.1, n))
    elseif field == :Ms
        uniform ? set_Ms(sim, 8e5) : set_Ms(sim, fill(8e5, n))
        add_exch(sim, 1.3e-11)
        add_anis(sim, 1e5)
        add_zeeman(sim, (0, 0, 1e5))
    elseif field == :Ku
        set_Ms(sim, 8e5)
        add_exch(sim, 1.3e-11)
        uniform ? add_anis(sim, 1e5) : add_anis(sim, fill(1e5, n))
        add_zeeman(sim, (0, 0, 1e5))
    else # :D
        set_Ms(sim, 8e5)
        add_exch(sim, 1.3e-11)
        uniform ? add_dmi(sim, 1e-3; type="interfacial") :
                  add_dmi(sim, fill(1e-3, n); type="interfacial")
        add_zeeman(sim, (0, 0, 1e5))
    end

    init_m0(sim, (1, 0.1, 0.05))
    run_sim(sim; steps=20, dt=1e-12, save_data=false, save_m_every=-1)
    return Array(sim.spin)
end

function relax_case(uniform::Bool)
    mesh = FDMesh(nx=16, ny=16, nz=2, dx=3e-9, dy=3e-9, dz=3e-9)
    sim = Sim(mesh; name="equiv_relax_$(uniform ? "fill" : "dense")", driver="LLG")
    set_Ms(sim, 8e5)
    add_exch(sim, 1.3e-11)
    add_zeeman(sim, (0, 0, 1e5))
    uniform ? set_alpha(sim, 0.1) : set_alpha(sim, fill(0.1, sim.n_total))
    init_m0(sim, (0.8, 0.6, 0.1))
    relax(sim; max_steps=60, stopping_dmdt=0.01, save_m_every=-1)
    return Array(sim.spin)
end

function test_fill_equivalence()
    set_precision(Float64)
    for field in (:alpha, :Ms, :Ku, :D)
        spin_fill = run_case(field, true)
        spin_dense = run_case(field, false)
        @test !iszero(maximum(abs.(spin_fill)))  # sanity: the system actually moved
        @test spin_fill == spin_dense            # bit-identical trajectories
    end

    m_fill = relax_case(true)
    m_dense = relax_case(false)
    @test m_fill == m_dense
end

@using_gpu()
test_functions("Fill equivalence", test_fill_equivalence, precisions=[Float64])
