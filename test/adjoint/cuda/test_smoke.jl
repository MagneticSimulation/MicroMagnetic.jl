#GPU smoke tests: prove the CUDA path the adjoint layer will stand on before any
#VJP kernel lands -- sim buffers on GPU, effective_field (exch/ani/demag/DMI)
#consistent with CPU, the cuFFT demag path, pbc2d demag and a relax end-to-end.
#Float64 is primary (the adjoint layer is Float64-only); one Float32 case for
#coverage of the C2R cuFFT branch.

#Switch the constructor-time backend without the "existing Sim(s)" warning noise
#(same runner trick as test/test_utils.jl).
function on_backend(f, name)
    saved_n_sims = MicroMagnetic._n_sims[]
    MicroMagnetic._n_sims[] = 0
    ok = MicroMagnetic.set_backend(name)
    MicroMagnetic._n_sims[] = saved_n_sims
    ok || error("backend $name is not available")
    return f()
end

function test_buffers_on_gpu()
    set_precision(Float64)
    sim = on_backend("cuda") do
        Sim(FDMesh(; nx=4, ny=4, nz=2, dx=5e-9, dy=5e-9, dz=5e-9); save_data=false)
    end
    #sim.spin/sim.field are created through create_zeros, which follows the
    #constructor-time default_backend.
    @test sim.spin isa CuArray{Float64}
    @test sim.field isa CuArray{Float64}
    @test sim.energy isa CuArray{Float64}
end

function build_stack_sim()
    mesh = FDMesh(; nx=6, ny=5, nz=2, dx=4e-9, dy=3e-9, dz=5e-9)
    sim = Sim(mesh; save_data=false)
    set_Ms(sim, 8e5)
    init_m0(sim, (0.3, -0.4, 0.5))
    add_exch(sim, 1.3e-11)
    add_anis(sim, 5e4; axis=(0, 0, 1))
    add_demag(sim)
    add_dmi(sim, 3e-3)  #bulk DMI: the one K^T != K term of the VJP table
    return sim
end

function test_effective_field_cross_backend()
    set_precision(Float64)
    cpu = on_backend("cpu") do
        s = build_stack_sim()
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    gpu = on_backend("cuda") do
        s = build_stack_sim()
        @test s.field isa CuArray
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    rel = maximum(abs.(cpu .- gpu)) / maximum(abs.(cpu))
    @test rel < 1e-9
end

#Regression anchor against the CPU/F incident values locked in test/test_demag.jl:
#same mesh/params, run on CUDA (cuFFT in-place path).
function test_demag_regression()
    set_precision(Float64)
    field = on_backend("cuda") do
        mesh = FDMesh(; nx=3, ny=2, nz=1)
        s = Sim(mesh; save_data=false)
        set_Ms(s, 8.6e5)
        init_m0(s, (0.1, 0.2, 1))
        add_demag(s)
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    expected = [-9615.99019074, -39898.78767025, -430282.70478141, -8664.33293854,
                -51323.59117349, -496012.77748287, -27749.99302205, -48965.78908591,
                -430282.70478141, -27749.99302205, -48965.78908591, -430282.70478141,
                -8664.33293854, -51323.59117349, -496012.77748287, -9615.99019074,
                -39898.78767025, -430282.70478141]
    @test isapprox(field, expected)
end

#xy-periodic demag: uniform in-plane M carries no bound charge (H ≈ 0), uniform
#open-axis M hits the DC column (H = -M along z, ≈ 0 elsewhere).  Anchors and
#tolerances mirror test/test_demag_pbc2d.jl's test_pbc2d_uniform (truncated
#image-sum residuals put the Float64 bound at 1e-3 * Ms, not machine precision).
function test_demag_pbc2d()
    set_precision(Float64)
    Ms = 8e5
    tol = 1e-3 * Ms
    nx, ny, nz = 4, 3, 2

    h_plane = on_backend("cuda") do
        mesh = FDMesh(; nx=nx, ny=ny, nz=nz, dx=2e-9, dy=3e-9, dz=4e-9, pbc="xy")
        s = Sim(mesh; save_data=false)
        set_Ms(s, Ms)
        init_m0(s, (0.6, -0.8, 0); norm=false)
        add_demag(s)
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    @test maximum(abs, h_plane) <= tol

    h_open = on_backend("cuda") do
        mesh = FDMesh(; nx=nx, ny=ny, nz=nz, dx=2e-9, dy=3e-9, dz=4e-9, pbc="xy")
        s = Sim(mesh; save_data=false)
        set_Ms(s, Ms)
        init_m0(s, (0, 0, 1); norm=false)
        add_demag(s)
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    target = zeros(3 * nx * ny * nz)
    target[3:3:end] .= -Ms
    @test maximum(abs, h_open .- target) <= tol
end

function test_relax_on_gpu()
    set_precision(Float64)
    m = on_backend("cuda") do
        mesh = FDMesh(; nx=4, ny=4, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
        sim = Sim(mesh; save_data=false)
        set_Ms(sim, 8e5)
        #off the hard-axis plane and demag-free: with m exactly perpendicular to
        #the easy axis the anisotropy torque vanishes identically and relax stalls
        #(and with a thin-film demag the equilibrium is in-plane, not ±z)
        init_m0(sim, (0.6, 0, -0.8))
        add_exch(sim, 1.3e-11)
        add_anis(sim, 5e4; axis=(0, 0, 1))
        relax(sim; max_steps=5000, stopping_dmdt=1e-5)
        Array(sim.spin)
    end
    mx, my, mz = m[1:3:end], m[2:3:end], m[3:3:end]
    @test all(abs.(mz) .> 1 - 1e-6)     #relaxed onto the ±z easy axis
    @test maximum(abs.(vcat(mx, my))) < 1e-4   #~stopping_dmdt=1e-5 residual tail
end

function test_float32()
    set_precision(Float32)
    cpu = on_backend("cpu") do
        s = build_stack_sim()
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    gpu = on_backend("cuda") do
        s = build_stack_sim()
        @test s.field isa CuArray{Float32}
        MicroMagnetic.effective_field(s, s.spin, 0.0)
        Array(s.field)
    end
    rel = maximum(abs.(cpu .- gpu)) / maximum(abs.(cpu))
    @test rel < 1e-4
end

test_buffers_on_gpu()
test_effective_field_cross_backend()
test_demag_regression()
test_demag_pbc2d()
test_relax_on_gpu()
test_float32()

#leave the constructor-time state as we found it (Float64 + CPU) in case this
#file is included from a larger suite.
set_precision(Float64)
on_backend(() -> nothing, "cpu")
