```@meta
ShareDefaultModule = true
```

# 标准问题 5（sim_with）

微磁标准问题 5 模拟矩形薄膜中的磁涡旋在自旋极化电流驱动下的动力学。
该问题用于检验微磁模拟捕捉涡旋结构复杂行为（包括运动与形变）的精度。本例中的矩形薄膜
尺寸为 100 nm × 100 nm × 10 nm，模拟追踪涡旋核位置随时间的演化——它由交换相互作用、
退磁场和自旋转移矩的共同效应驱动。

````@example
using MicroMagnetic
using Printf
using CairoMakie

@using_gpu()
````

定义一个粗略初始化涡旋的函数。

````@example
function init_fun(i, j, k, dx, dy, dz)
    x = i - 10
    y = j - 10
    r = sqrt(x^2 + y^2)
    if r < 2
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end
````

定义模拟参数。

````@example
args = (name="std5_sw", task_s=["relax", "dynamics"],      # 要执行的任务列表
        driver_s=["SD", "LLG"],                # 要使用的驱动列表
        mesh=FDMesh(; nx=20, ny=20, nz=2, dx=5e-9, dy=5e-9, dz=5e-9), # 网格配置
        Ms=8e5,                               # 饱和磁化强度
        A=1.3e-11,                            # 交换常数
        demag=true,                           # 启用退磁场
        m0=init_fun,                          # 初始磁化函数
        alpha=0.1,                            # Gilbert 阻尼参数
        stt=(model=:zhang_li, b=-72.438, J=(1, 0, 0), xi=0.05),  # 自旋转移矩
        steps=160,                            # 动力学步数
        dt=0.05ns,                            # 时间步长
        stopping_dmdt=0.01,                   # 弛豫停止判据
        saver_item=SaverItem(("Rx", "Ry"), ("<m>", "<m>"), compute_guiding_center),    #涡旋中心追踪
        dynamic_m_every=1);
nothing #hide
````

用 `sim_with` 函数运行模拟。

````@example
sim_with(args);
nothing #hide
````

生成涡旋动力学动画。使用 `LLG` 驱动时，ovf 文件存放在 `std5_sw_LLG` 文件夹：

````@example
ovf2movie("std5_sw_LLG"; output="../public/std5.mp4", component='x');
nothing #hide
````
![](../public/std5.mp4)
