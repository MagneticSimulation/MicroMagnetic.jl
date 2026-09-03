```@meta
ShareDefaultModule = true
```

# 用 NEB 模拟斯格明子塌缩

````@example
using CairoMakie
using DelimitedFiles
using CubicSplines
using MicroMagnetic
````

本例用 NEB 方法计算斯格明子塌缩到铁磁态的能量势垒。
首先用底层 API 描述所研究的体系。例如，体系是一个
带周期边界条件的薄膜（120x120x2 nm^3），考虑三种能量（交换、DMI 和 Zeeman）。

````@example
mesh = FDMesh(; nx=60, ny=60, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")

# 一个辅助函数，创建模拟并添加所有能量项。
function init_sim(; m0=nothing)
    sim = Sim(mesh; name="skx", driver="SD", save_data=false)
    set_Ms(sim, 3.84e5)
    add_exch(sim, 3.25e-12)
    add_dmi(sim, 5.83e-4)
    add_zeeman(sim, (0, 0, 120 * mT))
    if m0 !== nothing
        init_m0(sim, m0)
    end
    return sim
end
````

使用 NEB 可以分为两个阶段。第一阶段是准备初态和末态。我们假设
初态是一个磁性斯格明子，末态是铁磁态。

在这个方法中，我们将得到一个磁性斯格明子。斯格明子态保存为 'skx.vts'。

````@example
function relax_skx()
    function m0_fun_skx(i, j, k, dx, dy, dz)
        r2 = (i - 30)^2 + (j - 30)^2
        if r2 < 10^2
            return (0.01, 0, -1)
        end
        return (0, 0, 1)
    end

    sim = init_sim(; m0=m0_fun_skx)
    relax(sim; max_steps=2000, stopping_dmdt=0.01)
    save_vtk(sim, "skx")

    return sim
end
````

调用 relax_skx 方法得到磁性斯格明子态，并绘制磁化。

````@example 
sim = relax_skx()
plot_m(sim)
````

下面是第二阶段。

我们需要定义初态和末态，存放在 init\_images 列表中。
注意任何可接受的对象（如函数、元组或数组）都可以使用。
此外，如果有中间态，init\_images 列表也可以包含它。

```julia
init_images = [read_vtk("skx.vts"), (0, 0, 1)];
```

还需要一个插值数组来指定 NEB 模拟使用多少个 image。
注意插值数组的长度是 init\_images 的长度减一。例如，
若 init\_images = [read_vtk("skx.vts"), read_vtk("skx2.vts"), (0, 0, 1)]，则插值数组长度应为 2，
即类似 interpolation = [5,5]。

```julia 
interpolation = [6];
```

要使用 NEB，用同一个辅助函数创建一个 Sim 实例。

```julia
sim = init_sim();
```

创建 NEB 实例并设置 spring_constant，驱动可以是 "SD" 或 "LLG"

```julia 
neb = NEB(sim, init_images, interpolation; name="skx_fm", driver="SD");
````

neb.spring_constant = 1e7

弛豫整个体系

```julia
relax(neb; stopping_dmdt=0.1, save_vtk_every=1000, max_steps=5000)
```

运行模拟后，会生成能量文本文件（'skx_fm_energy.txt'）和对应的
距离文本文件（'skx_fm_distance.txt'）。

定义一个函数提取绘图所需的数据。

````@example 
function extract_data(; id=1)
    energy = readdlm("assets/skx_fm_energy.txt"; skipstart=2)
    dms = readdlm("assets/skx_fm_distance.txt"; skipstart=2)
    xs = zeros(length(dms[1, 1:end]))
    for i in 2:length(xs)
        xs[i] = sum(dms[id, 2:i])
    end

    et = energy[id, 2:end]
    e0 = minimum(et)
    energy_eV = (et .- e0) / meV

    spline = CubicSpline(xs, energy_eV)

    xs2 = range(xs[1], xs[end], 100)
    energy2 = spline[xs2]

    return xs, energy_eV, xs2, energy2
end

function plot_energy()
    fig = Figure(; resolution=(800, 480))
    ax = Axis(fig[1, 1]; xlabel="Distance (a.u.)", ylabel="Energy (meV)")

    xs, energy, xs2, energy2 = extract_data(; id=1)
    scatter!(ax, xs, energy; markersize=6, label="Initial energy path")
    lines!(ax, xs2, energy2)

    xs, energy, xs2, energy2 = extract_data(; id=500)
    scatter!(ax, xs, energy; markersize=6, label="Minimal energy path")
    lines!(ax, xs2, energy2)
    #linescatter!(ax, data[:,2]*1e9, data[:,5], markersize = 6)
    #linescatter!(ax, data[:,2]*1e9, data[:,6], markersize = 6)

    axislegend()

    save("energy.png", fig)

    return fig
end

plot_energy()
````

