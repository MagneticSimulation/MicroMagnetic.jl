```@meta
ShareDefaultModule = true
```

# Stoner–Wohlfarth 模型

链接：<https://en.wikipedia.org/wiki/Stoner%E2%80%93Wohlfarth_model>

当所研究体系的尺寸小于交换长度时，体系磁化是均匀的，可以用单个三维矢量描述。此时退磁场
可以用退磁张量与磁化简单相乘得到，退磁能等价于一个有效各向异性。
设外场 $\vec{H}=(H, 0, 0)$、易轴 $\hat{u}=(\cos\theta, \sin\theta, 0)$、单位磁化矢量
$\hat{m}= (\cos\phi, \sin\phi, 0)$，则有
```math
\begin{align}
E &= - K (\vec{m} \cdot \hat{u})^2 - \mu_0 M_s H \cos \phi \\
&= -\frac{K}{2}\left[ 1 + \cos (2(\theta-\phi))+ 4 h \cos\phi \right ]
\end{align}
```
其中 $h=H/H_k$，$H_k = 2K/(\mu_0 M_s)$。平衡态下能量对磁化方向的一阶导为零，即
```math
\frac{\partial E}{\partial \phi} = K [-\sin (2(\theta-\phi)) + 2 h \sin \phi] = 0
```
原则上每组给定的 $h$ 和 $\theta$ 都有对应的 $\phi$ 解。把 $\cos\phi$ 作为 $h$ 的函数作图即可得到
磁滞回线。翻转场可由能量对磁化方向的二阶导也为零额外求得，即
```math
\frac{\partial^2 E}{\partial \phi^2} = 2 K [\cos (2(\theta-\phi)) +  h \cos \phi] = 0
```
得到的翻转场为
```math
h_s=\frac{\left(1-t^2+t^4\right)^{1 / 2}}{1+t^2}
```
其中 $t=\tan ^{1 / 3} \theta$。特别地，当 $\theta=\pi/4$ 时 $h_s=1/2$。本例用 MicroMagnetic
验证这一结果。

我们取立方样品，因此退磁张量为 $N_x=N_y=N_z=1/3$，
即退磁本身对有效各向异性无贡献，所以本模拟忽略了退磁场。MicroMagnetic 脚本如下：

````@example
using MicroMagnetic
````

创建 4nm × 4nm × 4nm 立方几何的网格

````@example
mesh = FDMesh(; nx=4, ny=4, nz=4, dx=1e-9, dy=1e-9, dz=1e-9);
nothing #hide
````
创建使用 LLG 驱动和 BS23 积分器的模拟实例。
````@example
sim = Sim(mesh, driver="LLG", integrator="BS23", name="SW")

sim.driver.alpha = 0.5
sim.driver.integrator.tol = 1e-7

set_Ms(sim, 1e6)
init_m0(sim, (-1,1,0))
    
add_exch(sim, 1.3e-11)
add_anis(sim, 5e4, axis=(1, 1, 0))
add_zeeman(sim, (0,0,0))
nothing #hide
````

用 `hysteresis` 函数模拟磁滞回线。
````@example
Hs = [i*mT for i=-100:5:100]
hysteresis(sim, Hs, direction=(1,0,0), full_loop=false, stopping_dmdt=0.05, output="vts")
nothing #hide
````

对所用各向异性 $K_u=5\times 10^4$ A/m$^3$，预期翻转场为 $H_c = (1/2) H_K = 39788.7$ A/m。
用 `plot_ts` 函数绘制磁滞回线

````@example
using CairoMakie
fig = plot_ts("SW_llg.txt", x_key="Hx", ["m_x", "m_y"], x_unit=1/mT, xlabel="H (mT)", ylabel="m", mirror_loop=true);
````

```@setup
save("../public/sw.png", fig)
```
