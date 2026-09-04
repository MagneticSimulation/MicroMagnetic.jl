# 单位制

MicroMagnetic.jl 的模拟默认使用 **SI 单位制**。

## 微磁模型

在微磁模型中，所有参数均使用 SI 单位表示。下表列出了常用物理量及其对应单位。

| 物理量                       | 单位    | 用法示例                                                       |
|:---------------------------:|:-------:|:-------------------------------------------------------------:|
| **长度**                    | m       | `mesh = FDMesh(dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1)` |
| **时间**                    | s       | `run_sim(sim; steps=100, dt=1e-11)`                           |
| **饱和磁化强度**            | A/m     | `set_Ms(sim, 8.6e5)`                                          |
| **交换常数**                | J/m     | `add_exch(sim, 1.3e-11)`                                      |
| **DMI 常数**                | J/m²    | `add_dmi(sim, 1e-3)`                                          |
| **各向异性常数**            | J/m³    | `add_anis(sim, 1e5; axis=(0,0,1))`                            |
| **外加磁场**                | A/m     | `add_zeeman(sim, (-24.6mT, 4.3mT, 0))`                        |
| **旋磁比**                  | m/(A·s) | `sim.driver.gamma = 2.21e5`                                   |
| **温度**                    | K       | `add_thermal_noise(sim, 100.0)`                               |

## 原子模型

原子模型在适用的地方同样使用 SI 单位，但许多量通常以能量或磁矩的形式表示。

| 物理量                     | 单位      | 用法示例                                                       |
|:-------------------------:|:--------:|:-------------------------------------------------------------:|
| **长度**                  | m        | `mesh = CubicMesh(dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1)` |
| **时间**                  | s        | `run_sim(sim; steps=100, dt=1e-11)`                           |
| **磁矩**                  | A·m²     | `set_mu_s(sim, mu_s_1)`                                       |
| **交换相互作用**          | J        | `add_exch(sim, 50 * k_B)`                                     |
| **DMI 常数**              | J        | `add_dmi(sim, 5 * k_B)`                                       |
| **各向异性能量**          | J        | `add_anis(sim, 5 * k_B; axis=(0,0,1))`                        |
| **外加磁场**              | T        | `add_zeeman(sim, (0, 0, 0.1))`                                |
| **旋磁比**                | rad/(T·s) | `sim.driver.gamma = 1.76e11`                                  |
| **温度**                  | K        | `add_thermal_noise(sim, 100.0)`                               |

## 无量纲单位

在某些情况下，可以使用**无量纲单位**来简化模拟，这在理论研究或标度分析中特别有用。

### 无量纲化设置示例

在 $J = 1\,\mathrm{meV}$、$S = 1$、$a = 0.5\,\mathrm{nm}$ 的约定下,可用下表把
模拟结果换算回物理单位:

|物理量 | 换算系数 |  物理值  |
| :----:   | :----: | :----: |
| 距离 x   | $\hat{x}=a$ | 0.5nm |
| 时间 t   | $\hat{t}=\hbar S/J$ | 0.66 ps |
| 磁场 H | $\hat{H} = J/(\hbar \gamma S)$ | 8.63T |
| 速度 v | $\hat{v} = Ja/(\hbar S)$ | 759.63 m/s |
| 频率 $\omega$ | $\hat{\omega} = J/(\hbar S)$ | 1519.3 GHz |
| 温度 T | $\hat{T} = J/k_B$  |  11.6 K |

| 参数                    | 取值           |
|:-----------------------:|:--------------:|
| **晶格常数**            | $a = 1$        |
| **自旋长度**            | $S = 1$       |
| **旋磁比**              | $ \gamma = 1 $ |
| **磁矩**                | $ \mu_s = 1 $ |
| **交换常数**            | $ J = 1 $    |
| **DMI 强度**            | $ D/J = 0.09 $ |
| **磁场**                | $ H_0 = 0.00729 $ |
| **Gilbert 阻尼**        | $ \alpha = 0.04 $ |


下面是在 MicroMagnetic.jl 中使用无量纲单位搭建模拟的示例：

```julia
using MicroMagnetic

# 用归一化晶格单位定义立方网格
mesh = CubicMesh(; nx=200, ny=200, nz=1, dx=1, dy=1, dz=1)

# 用 Landau-Lifshitz-Gilbert (LLG) 驱动初始化模拟
sim = Sim(mesh; driver="LLG", name="test")
sim.driver.gamma = 1
sim.driver.alpha = 0.04

# 设置磁矩
set_mu_s(sim, 1)

# 初始化磁化方向
init_m0(sim, (1, 0.2, 0))

# 添加交换相互作用和 DMI
add_exch(sim, 1; name="exch")
add_dmi(sim, 0.09; name="dmi")

# 施加无量纲外磁场
add_zeeman(sim, (0, 0, 0.00729))
```
