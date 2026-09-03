```@meta
ShareDefaultModule = true
```

# 磁性斯格明子

本示例演示如何用 MicroMagnetic.jl 得到一个斯格明子。参数取自 PRL **111** 067203 (2013)。

!!! note "模拟所用参数"
    |参数 | 取值  |
    | :----:   | :----: |
    | 晶格常数 | $a = 0.5$ nm |
    | 自旋长度      | $S = 1$      |
    | 磁矩  |  $\mu_s = 2 \mu_B$ |
    | 交换常数 |  $J = 50 k_B$   |
    | DMI         | $D/J = 0.5$  |     |
    | 外磁场  | $H \mu_s /J  = 0.2$ |

定义一个函数来描述问题。

````@example
using MicroMagnetic
using CairoMakie

function m0_fun(i, j, k, dx, dy, dz)
    r2 = (i - 25)^2 + (j - 25)^2
    if r2 < 10^2
        return (0.01, 0, -1)
    end
    return (0, 0, 1)
end


function relax_system()
    mesh = CubicMesh(; nx=50, ny=50, nz=1, pbc="xy")

    #用 'SD' 驱动创建模拟
    sim = Sim(mesh; driver="SD", name="skx")

    set_mu_s(sim, mu_s_1) # 设置体系的 mu_s

    #用 `m0_fun` 函数初始化体系
    init_m0(sim, m0_fun)

    J = 50 * k_B
    add_exch(sim, J; name="exch")
    add_dmi(sim, 0.5 * J; name="dmi")

    Hz = 0.2 * J / mu_s_1
    add_zeeman(sim, (0, 0, Hz)) # Hz 的单位是特斯拉

    #弛豫体系
    relax(sim; max_steps=2000, stopping_dmdt=0.01)

    #把磁化保存到 vtk 文件
    save_vtk(sim, "skx"; fields=["exch", "dmi"])

    return sim
end

sim = relax_system();
nothing #hide
````

得到斯格明子后，绘制它

````@example
plot_m(sim)
````

