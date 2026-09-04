```@meta
ShareDefaultModule = true
```

# 动力学磁化率

动力学磁化率 $\chi$，也称 AC 磁化率，定义为样品磁化对磁场小变化的响应：
```math
\chi = \frac{d M}{d H}
```

施加交变磁场
```math
H(t) = H_0 + h \cos(\omega t)
```

$H(t)$ 引起的磁化为
```math
M(t) = M_0 + m \cos(\omega t - \phi)
```
其中 $\phi$ 是描述 $M(t)$ 相对 $H(t)$ 滞后的相位角。于是
```math
M(t) = M_0 + \chi' h \cos(\omega t) +  \chi'' h \sin(\omega t)
```
其中 $\chi'=\frac{m}{h}\cos\phi$，$\chi''=\frac{m}{h}\sin\phi$。可以用傅里叶变换计算 $\chi$：

```math
\chi(\omega) = \hat{m}(\omega)/ \hat{h}(\omega)
```
其中 $\hat{m}(\omega)$ 和 $\hat{h}(\omega)$ 分别是 $m(t)$ 和 $h(t)$ 的傅里叶变换。

````@example
using MicroMagnetic
using DelimitedFiles
using CairoMakie
using FFTW
````

启用 GPU 加速

````@example
@using_gpu()

mesh = FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);
nothing #hide
````

## 第一步：弛豫体系

模拟的第一步是弛豫体系，得到通常称为 "S" 态的稳定磁化组态。我们把这一过程封装在 `relax_system` 函数中。

````@example
function relax_system(mesh)
   sim = Sim(mesh; driver="SD", name="std4_bar")

   Ms = 8e5
   A = 1.3e-11

   set_Ms(sim, Ms)   # 设置饱和磁化强度
   add_exch(sim, A)  # 添加交换相互作用
   add_demag(sim)    # 添加退磁场

   init_m0(sim, (1, 0.25, 0.1))  # 初始化磁化
   relax(sim; stopping_dmdt=0.001)  # 弛豫体系

   return sim
end

sim = relax_system(mesh);
plot_m(sim; component='x')
````

## 第二步：施加外场

````@example
function time_fun(t)
    w = 2*pi*2.0e9
    return sinc(w*t)
end

function run_dynamics(sim)
    set_driver(sim; driver="LLG", alpha=0.015, gamma=2.211e5)
    sim.driver.integrator.tol = 1e-8
    add_zeeman(sim, (0, 500, 0), time_fun, name="zee")  # 沿 y 方向施加外磁场

    run_sim(sim; steps=10000, dt=1e-12, save_m_every=-1)  # 运行 10000 步模拟
end

function compute_chi(Ms=8e5)
    data, units = read_table("std4_bar_llg.txt")

    time = data["time"]
    dt = time[2] - time[1]
    N = length(time)
    println(dt)

    freq = fftshift(fftfreq(N, 1/dt))

    M = data["m_y"]
    H = data["zee_Hy"]
    fH = fftshift(fft(H))
    fM = fftshift(fft(M .- M[1]))

    a = real(fH)
    b = imag(fH)
    c = real(fM)
    d = imag(fM)

    rx = (a .* c .+ b .* d) ./ (a .* a .+ b .* b)
    ix = (b .* c .- a .* d) ./ (a .* a .+ b .* b)

    return freq*1e-9, rx*Ms, ix*Ms
end

function plot_chi()

    data = readdlm("assets/chi.txt")

    freq = data[:, 1]
    ix = data[:, 3]

    fig = Figure(; size=(400, 280), backgroundcolor=:transparent)
    ax = Axis(fig[1, 1]; xlabel="Frequency (GHz)", ylabel=L"\chi_y", backgroundcolor=:transparent)

    scatterlines!(ax, freq, ix; markersize=6, color=:blue, markercolor=:orange)
    xlims!(ax, 0, 15)
    ylims!(ax, -10, 200)
    return fig
end

if !isfile("assets/chi.txt")

    run_dynamics(sim)

    freq, rx, ix = compute_chi()

    s = div(length(freq), 2)
    data = [freq[s:end] rx[s:end] ix[s:end]]
    writedlm("assets/chi.txt", data)
end

fig = plot_chi()
````

```@setup
save("../public/chi.png", fig)
```

