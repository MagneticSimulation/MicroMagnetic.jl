```@meta
ShareDefaultModule = true
```

# 标准问题 5

链接：<https://www.ctcms.nist.gov/~rdm/std5/spec5.xhtml>


本脚本演示用 MicroMagnetic.jl 模拟矩形磁性薄膜，薄膜尺寸为 100 nm × 100 nm × 10 nm。
模拟分为两个主要步骤：
1. 弛豫体系，得到能量最优的磁化组态。
2. 模拟涡旋在电流驱动下的动力学。

````@example
using MicroMagnetic
using CairoMakie
using Printf
using DelimitedFiles

@using_gpu()
````

体系为 100 nm × 100 nm × 10 nm 的矩形磁性薄膜。

````@example
mesh = FDMesh(; nx=20, ny=20, nz=2, dx=5e-9, dy=5e-9, dz=5e-9);
nothing #hide
````

粗略初始化一个涡旋。

````@example
function init_fun(i, j, k, dx, dy, dz)
    x = i - 10
    y = j - 10
    r = (x^2 + y^2)^0.5
    if r < 2
        return (0, 0, 1)
    end
    return (-y / r, x / r, 0)
end
````

## 第一步：弛豫体系

该函数用最速下降驱动把体系弛豫到能量最优态。

````@example
function relax_system(mesh)
    sim = Sim(mesh; driver="SD", name="std5")

    A = 1.3e-11  # 交换常数
    Ms = 8e5     # 饱和磁化强度
    set_Ms(sim, Ms)
    add_exch(sim, A)  # 交换长度=5.7nm
    add_demag(sim)

    init_m0(sim, init_fun)
    relax(sim; max_steps=10000, save_m_every=-1)

    return sim
end

sim = relax_system(mesh);
nothing #hide
````

用 plot_m 函数绘制磁化分布。

````@example
plot_m(sim; arrows=(30, 30), component='x')
````

## 第二步：涡旋动力学

把驱动切回 LLG，并给体系加上自旋转移矩。

````@example
set_driver(sim; driver="LLG", alpha=0.1, gamma=2.211e5)
add_stt(sim, model=:zhang_li, P=1.0, Ms=8e5, xi=0.05, J=(1e12, 0, 0))

# 添加一个 SaverItem，每步保存涡旋导引中心
center = SaverItem(("cx", "cy"), ("<m>", "<m>"), compute_guiding_center)
run_sim(sim, steps=100, dt=5e-11, saver_item=center, save_m_every=1)
````

绘制涡旋中心随时间的变化。

````@example
fig = plot_ts("std5_llg.txt", ["cx", "cy"]; xlabel="Time (ns)", ylabel="Vortex center (nm)", y_units=[1e9, 1e9], transparency=true)
````

```@setup
save("../public/std5_center.png", fig)  #保存图像
```

````@example
ovf2movie("std5_LLG"; output="../public/std5.mp4", component='x');
nothing #hide
````
![](../public/std5.mp4)
