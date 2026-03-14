using MicroMagnetic
using Test
using NPZ

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts::Array)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function init_m0_function(i, j, k, dx, dy, dz)
    mx = sin(i/5)
    my = cos(i/5)
    mz = 0
    return (mx, my, mz)
end

function test_Laplaian()
    mesh = FDMesh(nx=10, ny=2, nz=2, dx=2e-9, dy=2e-9, dz=2e-9)
    sim = Sim(mesh; name="test", driver="SD")
    set_Ms(sim, 8e5)

    init_m0(sim, init_m0_function)

    exch = add_exch(sim, 1.3e-11)
    M = MicroMagnetic.build_exch_matrix(exch, sim)

    m = Array(sim.spin)
    m = reshape(m, (3, 40))
    mx = m[1, :]
    my = m[2, :]
    mz = m[3, :]

    fx = M*mx
    fy = M*my
    fz = M*mz
    MicroMagnetic.effective_field(sim, sim.spin)
    exch = Array(exch.field)
    f = reshape(exch, (3, 40))
    
    @test isapprox(fx , f[1, :], atol=1e-12)
    @test isapprox(fy , f[2, :], atol=1e-12)
    @test isapprox(fz , f[3, :], atol=1e-12)
end


function init_dw(i, j, k, dx, dy, dz)
    if i < 240
        return (1, 0.1, 0)
    elseif i < 260
        return (0, 1, 0.1)
    else
        return (-1, 0.1, 0)
    end
end

function relax_system()
    mesh = FDMesh(; nx=500, ny=1, nz=11, dx=2e-9, dy=2e-9, dz=1e-9)
    sim = Sim(mesh; name="relax_llg", driver="LLG", integrator="GPSM")
    sim.driver.precession = false
    sim.driver.alpha = 0.5
    
    set_Ms(sim, 8.6e5)

    add_exch(sim, 1.3e-11)
    add_anis(sim, 1e5; axis=(1, 0, 0))

    init_m0(sim, init_dw)
    relax(sim; max_steps=2000, stopping_dmdt=0.1)
    m = MicroMagnetic.average_m(sim)
    @test abs(m[1]) < 0.01
    @test abs(m[2]) < 0.04
    @test abs(m[3]) < 0.01
end


relax_system()

