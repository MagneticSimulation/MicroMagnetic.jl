```@meta
ShareDefaultModule = true
```

# 斯格明子晶格

本示例演示如何用 MicroMagnetic.jl 得到斯格明子晶格。参数取自 PRL **108** 017206 (2012)，为无量纲量。

````@example
using MicroMagnetic
using CairoMakie

function m0_fun(i, j, k, dx, dy, dz)
    i0, j0, r = 166, 96, 25
    i1 = i % i0
    j1 = j % j0

    if ((i1 - r)^2 + (j1 - r)^2 < r^2)
        return (0.05, 0.01, -1)
    elseif ((i1 - i0 / 2.0 - r)^2 + (j1 - j0 / 2.0 - r)^2 < r^2)
        return (0.05, 0.01, -1)
    end

    return (0, 0, 1)
end

function relax_system()
    mesh = CubicMesh(; nx=166 * 2, ny=96 * 3, nz=1, pbc="xy")
    sim = Sim(mesh; driver="SD", name="skx_latttice")
    set_mu_s(sim, 1.0)

    #设置交换、DMI 和 Zeeman
    add_exch(sim, 1.0; name="exch")
    add_zeeman(sim, (0, 0, 3.75e-3))
    add_dmi(sim, 0.09; name="dmi")

    init_m0(sim, m0_fun)
    relax(sim; max_steps=2000, stopping_dmdt=1e-5)

    save_vtk(sim, "skx_latttice.vts")

    return sim
end
````

调用 `relax_system` 函数得到斯格明子晶格。

````@example
sim = relax_system();
nothing #hide
````

得到斯格明子后，用以下脚本绘制

````@example
fig = plot_m(sim)
````
