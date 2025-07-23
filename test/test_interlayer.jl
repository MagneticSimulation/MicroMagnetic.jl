using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function m0_fun(i, j, k, dx, dy, dz)
    if k == 1
        return (1, 2, 3)
    end
    if k == 3
        return (4, 5, 6)
    end
    return (1, 1, 1)
end

function init_single_skx(i, j, k, dx, dy, dz)
    r = 20
    if k == 3
        if (i - 10)^2 + (j - 30)^2 < r^2
            return (0, -0.01, -1)
        end
        return (0, 0.01, 1)
    end
    if k == 1
        if (i - 10)^2 + (j - 30)^2 < r^2
            return (0, 0.01, 1)
        end
        return (0, -0.01, -1)
    end
end

function spatial_Ms(i, j, k, dx, dy, dz)
    if k == 1 || k == 3
        return 8e5
    end
    return 0
end

function test_interlayer()
    J = 1e-5
    Ms = 8e5
    dz = 2e-9
    Dx, Dy, Dz = 2e-3, 1e-4, 3e-5

    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=3, ny=5, nz=4, pbc="xy")
    sim = Sim(mesh)
    set_Ms(sim, spatial_Ms)
    init_m0(sim, m0_fun; norm=false)
    exch = add_exch_int(sim, J; k1=1, k2=3)
    dmi = add_dmi_int(sim, (Dx, Dy, Dz); k1=1, k2=3)

    MicroMagnetic.effective_field(sim, sim.spin)

    expected_exch = [4, 5, 6, 0, 0, 0, 1, 2, 3, 0, 0, 0] .* J / (mu_0 * Ms * dz)
    expected_exch = vcat([repeat(expected_exch[i:(i + 2)], 15)
                          for i in 1:3:length(expected_exch)]...)

    fx1 = 1 / (mu_0 * Ms * dz) * MicroMagnetic.cross_x(Dx, Dy, Dz, 4.0, 5.0, 6.0)
    fy1 = 1 / (mu_0 * Ms * dz) * MicroMagnetic.cross_y(Dx, Dy, Dz, 4.0, 5.0, 6.0)
    fz1 = 1 / (mu_0 * Ms * dz) * MicroMagnetic.cross_z(Dx, Dy, Dz, 4.0, 5.0, 6.0)
    fx2 = 1 / (mu_0 * Ms * dz) * MicroMagnetic.cross_x(1.0, 2.0, 3.0, Dx, Dy, Dz)
    fy2 = 1 / (mu_0 * Ms * dz) * MicroMagnetic.cross_y(1.0, 2.0, 3.0, Dx, Dy, Dz)
    fz2 = 1 / (mu_0 * Ms * dz) * MicroMagnetic.cross_z(1.0, 2.0, 3.0, Dx, Dy, Dz)

    expected_dmi = [fx1, fy1, fz1, 0, 0, 0, fx2, fy2, fz2, 0, 0, 0]
    expected_dmi = vcat([repeat(expected_dmi[i:(i + 2)], 15)
                         for i in 1:3:length(expected_dmi)]...)

    @test isapprox(expected_exch, Array(exch.field))
    @test isapprox(expected_dmi, Array(dmi.field))
end

function test_interlayer_exch()

    J = 1e-5

    function spatial_Js(i,j,k,dx,dy,dz)
        return J
    end
    Ms = 8e5
    dz = 2e-9

    nx, ny, nz = 30, 80, 3
    mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=nx, ny=ny, nz=nz, pbc="xy")
    sim = Sim(mesh)
    set_Ms(sim, spatial_Ms)
    init_m0(sim, init_single_skx)

    exch = add_exch_int(sim, spatial_Js; k1=1, k2=3)

    MicroMagnetic.effective_field(sim, sim.spin)

    f = Array(exch.field)
    f = reshape(f, (3, nx, ny, nz))
    f1 = f[:, :, :, 1]
    f3 = f[:, :, :, 3]

    s = f1 .+ f3

    @test maximum(s) == 0
    @test minimum(s) == 0
end

@using_gpu()

test_functions("Interlayer", test_interlayer, test_interlayer_exch)
