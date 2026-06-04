using MicroMagnetic

function vcross(u, v)
    return (u[2] * v[3] - u[3] * v[2],
            u[3] * v[1] - u[1] * v[3],
            u[1] * v[2] - u[2] * v[1])
end

function vscale(a, u)
    return (a * u[1], a * u[2], a * u[3])
end

function vnorm(u)
    n = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
    return (u[1] / n, u[2] / n, u[3] / n)
end

function setup_twin_monoclinic_example()
    mesh = FDMesh(; nx=160, ny=120, nz=6, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh; driver="SD", name="twin_monoclinic_anisotropy")

    Ms = 4.8e5
    A = 1.2e-11
    add_exch(sim, A)
    add_demag(sim)
    add_zeeman(sim, (0, 0, 1.0e4))

    c_axis = (0.0, 0.0, 1.0)
    b_left = vnorm((-1.0, 1.0, 0.0))
    b_right = vnorm((1.0, 1.0, 0.0))
    a_left = vnorm(vcross(b_left, c_axis))
    a_right = vnorm(vcross(b_right, c_axis))

    function axis_b_fun(x, y, z)
        return x < 400e-9 ? b_left : b_right
    end

    function axis_a_fun(x, y, z)
        return x < 400e-9 ? a_left : a_right
    end

    function axis_u111_fun(x, y, z)
        a = x < 400e-9 ? a_left : a_right
        return vnorm(vscale(sqrt(2 / 3), a) .+ vscale(1 / sqrt(3), c_axis))
    end

    set_Ms(sim, Ms)
    add_twin_mono_anis(sim; Ka=0.0, Kb=2.0e4, Kaa=0.0, Kbb=1.3e4, Kab=0.0,
                       Ku=8.0e3, axis_a=axis_a_fun, axis_b=axis_b_fun,
                       axis_u111=axis_u111_fun)

    function init_m(x, y, z)
        phase = sin(2pi * y / 180e-9)
        return phase > 0 ? (0.15, 0.1, 0.98) : (-0.15, -0.1, -0.98)
    end

    init_m0(sim, init_m)
    relax(sim; max_steps=5000, stopping_dmdt=0.02)
    save_vtk(sim, "twin_monoclinic_anisotropy")
    return sim
end

setup_twin_monoclinic_example()
