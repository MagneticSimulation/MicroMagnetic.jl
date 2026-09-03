# 方程

## 原子自旋模型

原子自旋模型的基本假设是每个格点关联一个磁矩 $\mu_s$。
对于轨道磁矩被猝灭的金属体系，磁矩主要与其自旋角动量相关：
```math
\mathbf{\mu} = - g \mu_B \mathbf{S} = - \hbar  \gamma \mathbf{S}
```
其中 $\mu_B=e \hbar /(2m)$ 是玻尔磁子，$e(>0)$ 是元电荷，
$\gamma=g\mu_B/\hbar (>0)$ 是旋磁比，$g=2$ 是 g 因子。
磁矩的动力学由 LLG 方程描述：

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H}_\mathrm{eff} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
```

其中 $\mathbf{m}$ 是磁矩的单位矢量，$\mathbf{H}_\mathrm{eff}$ 是有效场。

与连续的微磁模型不同，原子自旋模型是离散的，因此有效场 $\mathbf{H}_\mathrm{eff}$ 定义为总哈密顿量对 $\mathbf{m}$ 的偏导数，即

```math
\mathbf{H}_{\mathrm{eff}}=-\frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \mathbf{m}}
```

其中 $\mathcal{H}$ 是总哈密顿量，包括交换相互作用、Dzyaloshinskii-Moriya 相互作用、偶极相互作用、各向异性、Zeeman 相互作用等。


交换相互作用为

```math
\mathcal{H}_\mathrm{ex} = -J \sum_{\langle i, j\rangle} \mathbf{m}_{i} \cdot \mathbf{m}_{j}
```
其中 $J$ 是交换常数，$\langle i,j \rangle$ 表示格点 $i$ 与 $j$ 之间的一对唯一近邻，并且假设求和时每对只取一次。$J$ 为正时体系倾向于铁磁态，$J$ 为负时倾向于反铁磁态，这可以通过极小化交换哈密顿量看出。要数值求解 LLG 方程，需要每个格点上的有效场：

```math
\mathbf{H}_\mathrm{ex, i} = \frac{J}{\mu_s} \sum_{\langle i, j\rangle} \mathbf{m}_{j}.
```

类似地，Dzyaloshinskii-Moriya 相互作用为

```math
\mathcal{H}_\mathrm{dmi} =  \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right)
```
其中 $\mathbf{D}_{i j}$ 是 DM 矢量。典型情形有两种：(1) $\mathbf{D}_{i j} = D \hat{r}_{ij}$，对应体 DMI（Bulk DMI）；(2) $\mathbf{D}_{i j} = D \hat{r}_{ij} \times \hat{z}$，对应界面 DMI（interfacial DMI）。与交换相互作用不同，当两个相邻磁矩互相垂直时 DM 哈密顿量取极小。对应的有效场为
```math
\mathbf{H}_\mathrm{dmi, i} = \frac{1}{\mu_s} \sum_{\langle i, j\rangle} \mathbf{D}_{i j} \times \mathbf{m}_{j}.
```

另外两种常见的相互作用是各向异性和 Zeeman 项：

```math
\mathcal{H}_\mathrm{an} = - K \sum_{i}\left(\mathbf{m}_{i} \cdot \hat{u}\right)^{2} \\
\mathcal{H}_\mathrm{ze} =  - \mu_s \mathbf{m}_i \cdot \mathbf{H}
```

各向异性的有效场为
```math
\mathbf{H}_\mathrm{an, i} = \frac{2K}{\mu_s}  (\mathbf{m}_{i} \cdot \hat{u}) \hat{u}.
```

### 多铁绝缘体

对于 Cu$_2$OSeO$_3$ 等多铁材料，非共线自旋结构通过 $d-p$ 杂化机制诱导出局域电极化。局域电偶极矩 $\mathbf{P}$ 依赖于外加静态磁场相对于晶轴的方向。[PRL 128, 037201 (2022)]


##### $H_0//001$：
```math
\mathbf{P}_i=\lambda\left(-m_{i, z} m_{i, x}, m_{i, y} m_{i, z}, \frac{-m_{i, x}^2+m_{i, y}^2}{2}\right)
```

##### $H_0//110$：
```math
\mathbf{P}_i=\lambda\left(-m_{i, x} m_{i, y}, \frac{-m_{i, x}^2+m_{i, z}^2}{2}, m_{i, y} m_{i, z}\right)
```

##### $H_0//111$：

```math
\mathbf{P}_{i,x}=-\lambda \frac{m_{i, x}\left(\sqrt{2} m_{i, y}+m_{i, z}\right)}{\sqrt{3}} \\

\mathbf{P}_{i,y} = \lambda \frac{-m_{i, x}^2+m_{i, y}\left(m_{i, y}-\sqrt{2} m_{i, z}\right)}{\sqrt{6}} \\

\mathbf{P}_{i,z} =
-\lambda \frac{m_{i, x}^2+m_{i, y}^2-2 m_{i, z}^2}{2 \sqrt{3}}
```

与高频激光相关的哈密顿量为

```math
\mathcal{H}_\mathrm{laser} =  -\sum_{i} \mu_s \mathbf{m}_i \cdot \mathbf{H}(t) - \sum_{i} \mathbf{P}_i \cdot \mathbf{E}(t)
```

与电场相对应的有效场为

##### $H_0//001$：
```math
H_x = \frac{\lambda}{\mu_s}  (-E_z m_x  - E_x m_z )\\
H_y = \frac{\lambda}{\mu_s}  (E_z m_y  + E_y m_z )\\
H_z = \frac{\lambda}{\mu_s}  (-E_x m_x  + E_y m_y )\\
```

##### $H_0//110$：
```math
H_x = \frac{\lambda}{\mu_s}  (-E_y m_x  - E_x m_y )\\
H_y = \frac{\lambda}{\mu_s}  (-E_x m_x  + E_z m_z )\\
H_z = \frac{\lambda}{\mu_s}  (E_z m_y  + E_y m_z )\\
```

##### $H_0//111$：
```math
H_x = -\frac{\lambda}{\sqrt{3}\mu_s}  (\sqrt{2}(E_y m_x + E_x m_y)  +E_z m_x + E_x m_z )\\
H_y = -\frac{\lambda}{\sqrt{3}\mu_s}  (\sqrt{2}(E_x m_x - E_y m_y)  +E_z m_y + E_y m_z )\\
H_z = -\frac{\lambda}{\sqrt{3}\mu_s}  (E_x m_x  + E_y m_y - 2 E_z m_z)\\
```

## 微磁模型

在微磁学中，有效场可以由总微磁能量计算：

```math
\mathbf{H}_{\mathrm{eff}}=-\frac{1}{\mu_{0} M_{s}} \frac{\delta E}{\delta \mathbf{m}}
```

典型的能量项有

- **交换能**

```math
  E_\mathrm{ex} = \int_{V} A (\nabla \mathbf{m})^2 \mathrm{d}V
```

  其中 $(\nabla \mathbf{m})^{2}=\left(\nabla m_{x}\right)^{2}+\left(\nabla m_{y}\right)^{2}+\left(\nabla m_{z}\right)^{2}$。对应的有效场为

```math
  \mathbf{H}_{\mathrm{ex}}=\frac{2 A}{\mu_{0} M_{s}} \nabla^{2} \mathbf{m}
```

- **Zeeman 能**

```math
  E_\mathrm{ex} = -  \mu_0 \int_{V}  \mathbf{H} \cdot \mathbf{M} \mathrm{d}V
```

  如预期的那样，有效场就是 $\mathbf{H}$。

- **各向异性**

  单轴各向异性能为

```math
  E_\mathrm{anis} = -\int_{V} K_{u} (\mathbf{m} \cdot \hat{u})^2 \, dV
```

  由此可以计算有效场：

```math
  \mathbf{H}_{\mathrm{an}}=\frac{2 K_u}{\mu_0 M_s}\left(\mathbf{m} \cdot \hat{u}\right) \hat{u}
```

- **立方各向异性**

  立方各向异性能为

```math
  E_\mathrm{cubic} = -\int_{V} K_c (m_x^4 + m_y^4 + m_z^4) \, dV
```

  因此对应的有效场为

```math
  \mathbf{H}_{\mathrm{cubic}}= \frac{4 K_c}{\mu_0 M_s}  
  \left( m_x^3 \mathbf{e}_x + m_y^3 \mathbf{e}_y + m_z^3 \mathbf{e}_z \right)
```


- **六角各向异性**

  六角各向异性的能量密度为：

```math
E = K_1 \sin^2 \theta + K_2 \sin^4 \theta + K_3 \sin^6 \theta \cos 6\phi
```

这里 $\theta$ 是磁化矢量 $\mathbf{m}$ 与 $c$ 轴（$z$ 轴）的夹角。
$\phi$ 是 $\mathbf{m}$ 在六角面上的投影相对于 $x$ 轴的夹角。

利用恒等式 $\cos 6x = -\sin^6 x + 15 \cos^2 x \sin^4 x - 15 \cos^4 x \sin^2 x + \cos^6 x$，能量密度
可以用 $m_x$、$m_y$、$m_z$ 改写为：
```math
E = K_1 (1 - m_z^2) + K_2 (1 - m_z^2)^2 + K_3 \left( m_x^6 - 15m_x^4 m_y^2 + 15m_x^2 m_y^4 - m_y^6 \right)
```
因此，对应的有效场为：
```math
\mathbf{H}_\mathrm{eff} = -\frac{6K_3}{\mu_0 M_s} \left( m_x^5 - 10m_x^3m_y^2 + 5m_xm_y^4 \right) \mathbf{e}_x 
 -\frac{6K_3}{\mu_0 M_s} \left( -5m_x^4m_y + 10m_x^2m_y^3 - m_y^5 \right) \mathbf{e}_y 
+ \frac{2m_z}{\mu_0 M_s} \left[ K_1 + 2K_2(1 - m_z^2) \right] \mathbf{e}_z
```


### Dzyaloshinskii-Moriya 能量

在连续极限下，DMI 能量密度 $w_\mathrm{dmi}$ 与所谓的 *Lifshitz 不变量* 相关，其形式为

```math
L^{(k)}_{ij} = m_i \frac{\partial m_j}{\partial x_k} - m_j \frac{\partial m_i}{\partial x_k}.
```

DMI 能量密度的形式取决于对称性类。对于对应于对称性类 $T$ 或 $O$ 的体 DMI，表达式为
```math
  w_\mathrm{dmi} = D(L^{(z)}_{yx} + L^{(y)}_{xz} + L^{(x)}_{zy}) = D \mathbf{m} \cdot (\nabla \times \mathbf{m}).
```
相应的有效场为
```math
  \mathbf{H}_\mathrm{dmi}=-\frac{2D}{\mu_0 M_s} (\nabla \times \mathbf{m}).
```
对于具有界面 DMI 的薄膜或对称性类为 $C_{nv}$ 的晶体，能量密度为
```math
  w_\mathrm{dmi}=D (L_{x z}^{(x)}+L_{y z}^{(y)} )=D\left(\mathbf{m} \cdot \boldsymbol{\nabla} m_z-m_z \boldsymbol{\nabla} \cdot \mathbf{m}\right),
```
有效场为
```math
  \mathbf{H}_\mathrm{dmi}=-\frac{2 D}{\mu_0 M_s} (\mathbf{e}_y \times \frac{\partial \mathbf{m}}{\partial x} - \mathbf{e}_x \times \frac{\partial \mathbf{m}}{\partial y}).
```
对于对称性类为 $D_{2d}$ 的晶体，DMI 能量密度为 $w_{\mathrm{dmi}}=D (L_{x z}^{(y)}+L_{y z}^{(x)})$，由此得到有效场
```math
  \mathbf{H}_\mathrm{dmi}=-\frac{2 D}{\mu_0 M_s} (\mathbf{e}_y \times \frac{\partial \mathbf{m}}{\partial y} - \mathbf{e}_x \times \frac{\partial \mathbf{m}}{\partial x} ).
```
尽管不同对称性的有效场形式不同，数值实现可以统一写为
```math
  \mathbf{H}_\mathrm{dmi, i} = -\frac{1}{\mu_0 M_s} \sum_{j \in N_i} D_{ij} \frac{\mathbf{e}_{ij} \times \mathbf{m}_j}{\Delta_{i j}},
```
其中 $D_{ij}$ 是有效 DMI 常数，$\mathbf{e}_{ij}$ 是 DMI 矢量。
对于体 DMI，$\mathbf{e}_{ij} = \hat{\mathbf{r}}_{ij}$，其中 $\hat{\mathbf{r}}_{ij}$ 是单元 $i$ 与单元 $j$ 之间的单位矢量。
对于界面 DMI，$\mathbf{e}_{ij} = \mathbf{e}_z \times \mathbf{\hat{r}}_{ij}$，即 6 个近邻 $N_{i}=\{-x,+x,-y,+y,-z, +z\}$ 对应的
$\mathbf{e}_{ij}=\{-\mathbf{e}_y, \mathbf{e}_y, \mathbf{e}_x, -\mathbf{e}_x, 0, 0\}$。
对于对称性类 $D_{2d}$，$\mathbf{e}_{ij}=\{\mathbf{e}_x, -\mathbf{e}_x, -\mathbf{e}_y, \mathbf{e}_y, 0, 0\}$。
如果给出的是逐单元的 DMI，有效 DMI 常数可以计算为
```math
  D_{i j}=\frac{2 D_i D_j}{D_i+D_j}.
```

- **体 DMI 能量** 体 DMI 能量为

```math
  E_{\mathrm{dmi}} = \int_V D \mathbf{m} \cdot (\nabla \times \mathbf{m}) \, \mathrm{d}V
```

  因此有效场为

```math
  \mathbf{H}_\mathrm{D}=-\frac{2 D}{\mu_{0} M_{s}}(\nabla \times \mathbf{m})
```

- **静磁能**

```math
  E_{\mathrm{d}}=-\frac{\mu_{0}}{2} \int_{V} \mathbf{H}_{\mathrm{d}}(\mathbf{r}) \cdot \mathbf{M}(\mathbf{r}) d V
```

```math
  \mathbf{H}_{\mathrm{d}}(\mathbf{r})=\frac{1}{4 \pi}\left(\int_{V} \rho_{m}\left(\mathbf{r}^{\prime}\right) \frac{\mathbf{r}-\mathbf{r}^{\prime}}{\left|\mathbf{r}-\mathbf{r}^{\prime}\right|^{3}} \mathrm{d}^{3} r^{\prime}+\int_{S} \sigma_{m}\left(\mathbf{r}^{\prime}\right) \frac{\mathbf{r}-\mathbf{r}^{\prime}}{\left|\mathbf{r}-\mathbf{r}^{\prime}\right|^{3}} \mathrm{d}^{2} r^{\prime}\right)
```

## LLG 方程

`LLG` 驱动求解标准 LLG 方程，可写为

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
```

对应的 LL 形式为

```math
(1+\alpha^2)\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} - \alpha \gamma \mathbf{m} \times (\mathbf{m} \times \mathbf{H})
```

## 惯性 LLG 方程

惯性 LLG（Inertial LLG）方程描述带惯性效应的磁化动力学，它在经典 LLG 方程的基础上增加了二阶时间导数项：

```math
\frac{d \mathbf{m}}{dt} = - \gamma \mathbf{m}  \times \mathbf{H}_\mathrm{eff} +\alpha \mathbf{m} \times \frac{d \mathbf{m}}{dt} + \eta \mathbf{m} \times \frac{d^2 \mathbf{m}}{dt^2}
```

其中 $\mathbf{m}$ 是归一化磁化矢量（$|\mathbf{m}| = 1$），$\gamma$ 是旋磁比，$\mathbf{H}_\mathrm{eff}$ 是有效场，
$\alpha$ 是 Gilbert 阻尼参数，$\eta = \alpha \tau$ 是惯性参数，$\tau$ 是角动量弛豫时间。

为便于数值实现，引入速度变量：

```math
\mathbf{v} = \frac{d \mathbf{m}}{dt}
```

代入该定义后，ILLG 方程变为：

```math
\mathbf{v} = - \gamma \mathbf{m}  \times \mathbf{H}_\mathrm{eff} +\alpha \mathbf{m} \times \mathbf{v} + \eta \mathbf{m} \times \frac{d \mathbf{v}}{dt}
```

两边与 $\mathbf{m}$ 做叉积：

```math
\mathbf{m} \times \mathbf{v} =
-\gamma \mathbf{m} \times ( \mathbf{m} \times \mathbf{H}_\mathrm{eff})
+ \alpha \mathbf{m} \times (\mathbf{m} \times \mathbf{v})
+ \eta \mathbf{m} \times (\mathbf{m} \times \frac{d \mathbf{v}}{dt})
```

对惯性项应用矢量三重积恒等式 $\mathbf{a}\times (\mathbf{b} \times \mathbf{c}) = (\mathbf{a}\cdot\mathbf{c})\mathbf{b}-(\mathbf{a}\cdot\mathbf{b})\mathbf{c}$：

```math
\eta \mathbf{m} \times (\mathbf{m} \times \frac{d \mathbf{v}}{dt}) = 
\eta (\mathbf{m} \cdot \frac{d \mathbf{v}}{dt}) \mathbf{m}- \eta \frac{d \mathbf{v}}{dt}
= -\eta v^2 \mathbf{m} - \eta \frac{d \mathbf{v}}{dt}
```

化简 $(\mathbf{m} \cdot \frac{d \mathbf{v}}{dt}) = -v^2$ 来自对约束 $\mathbf{m} \cdot \mathbf{v} = 0$ 求导，而该约束本身来自磁化模长守恒 $|\mathbf{m}| = 1$。于是 ILLG 方程可以化为如下适于数值积分的一阶耦合方程组：

```math
\begin{aligned}
\frac{d \mathbf{m}}{dt} &= \mathbf{v} \\
\frac{d \mathbf{v}}{dt} &=
-\frac{1}{\eta}\mathbf{m} \times \mathbf{v} 
-\frac{\gamma}{\eta} \mathbf{m} \times ( \mathbf{m} \times \mathbf{H}_\mathrm{eff})
+ \frac{\alpha}{\eta}  \mathbf{m} \times (\mathbf{m} \times \mathbf{v}) - v^2 \mathbf{m}
\end{aligned}
```

## 自旋转移矩

自旋转移矩（spin transfer torque, STT）是自旋电子学中的基本现象：自旋极化的电流对磁矩施加力矩，从而实现对磁化状态的操控。该效应通常通过扩展 Landau-Lifshitz-Gilbert（LLG）方程来建模。

### Zhang-Li 模型

Zhang-Li 模型把 STT 纳入 LLG 方程。磁化矢量 **m** 的方程为：

```math
\frac{d \mathbf{m}}{dt} = -\gamma \mathbf{m} \times \mathbf{H}_\mathrm{eff} + \alpha \mathbf{m} \times \frac{d \mathbf{m}}{dt}
- b \mathbf{m} \times [\mathbf{m} \times (\mathbf{j} \cdot \nabla) \mathbf{m}] 
- \xi b \mathbf{m} \times (\mathbf{j} \cdot \nabla) \mathbf{m} 
```

其中 **j** 是电流密度矢量，$\xi$ 是非绝热参数。系数 $b$ 定义为：
```math
b = \frac{P \mu_B}{e M_s (1 + \xi^2)}
```
其中 $P$ 是自旋极化率，$\mu_B$ 是玻尔磁子，$e>0$ 是元电荷。

注意 $\mathbf{m} \times [ \mathbf{m} \times (\mathbf{j} \cdot \nabla) \mathbf{m}] = -(\mathbf{j} \cdot \nabla) \mathbf{m}$，方程可化简为：

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H}_\mathrm{eff}  + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
- (\mathbf{u} \cdot \nabla) \mathbf{m} + \beta [\mathbf{m}\times (\mathbf{u} \cdot \nabla)\mathbf{m}]
```

其中
```math
\mathbf{u} = -b \mathbf{J}, \, \, \, \, \beta = \xi
```
**u** 的单位是 m/s。

### 其他表述

其他模型可能用不同方式定义电流强度参数。作为比较，$\mathbf{u}$ 定义为：

```math
\mathbf{u} = -\frac{P g \mu_B}{2 e M_s} \mathbf{j}
```

关于模型差异的详细讨论见[该文献](https://www.ctcms.nist.gov/~rdm/std5/spec5.xhtml)。

### 原子模型

在原子模拟中，电流强度针对离散自旋做了适配：

```math
\mathbf{u} = -\frac{P g \mu_B a^3}{2 e \mu_s} \mathbf{j} = -\frac{P a^3}{2 e S} \mathbf{j}
```

其中 $a$ 是晶格常数，$\mu_s$ 是磁矩，$S = |\mathbf{S}|$ 是局域自旋矢量的模。

## Slonczewski 自旋转移矩

包含 Slonczewski 力矩的 LLG 方程为

```math
\frac{d\mathbf{m}}{dt} = -\gamma \mathbf{m} \times \mathbf{H}_{\text{eff}}
   + \alpha \left( \mathbf{m} \times \frac{d\mathbf{m}}{dt} \right)
   + \gamma \beta \epsilon \left[ \mathbf{m} \times (\mathbf{m}_p \times \mathbf{m}) \right]
   - \gamma \beta \xi \mathbf{m} \times \mathbf{m}_p
```

其中 $\mathbf{m}=\mathbf{M}/M_s$，$\gamma$ 是旋磁比，$\mathbf{m}_p$ 是电子极化方向，$\xi$ 是次级自旋力矩参数。
系数 $\beta$ 定义为：
```math
\beta =  \frac{\hbar}{\mu_0 e} \frac{J}{t M_s}
```
其中 $e$ 是元电荷（C），$J$ 是电流密度（A/m²），$t$ 是自由层厚度（m），$M_s$ 是饱和磁化强度（A/m）。系数 $\epsilon$ 定义为
```math
\epsilon = \frac{P \Lambda^2}{(\Lambda^2 + 1) + (\Lambda^2 - 1)(\mathbf{m} \cdot \mathbf{m}_p)}
```
其中 $P$ 是自旋极化率，$\Lambda$ 是 Slonczewski 参数。

当 $\epsilon=P/2$ 为常数（对应 $\Lambda=1$）时，Slonczewski 力矩退化为 [PRL **102**, 067206 (2009)]
```math
\frac{\partial \mathbf{m}}{\partial t} = -\gamma \, \mathbf{m} \times \mathbf{H}_{\text{eff}} 
+ \alpha \, \mathbf{m} \times \frac{\partial \mathbf{m}}{\partial t} 
- a_J \, \mathbf{m} \times (\mathbf{m} \times \mathbf{m}_p) 
- b_J \, \mathbf{m} \times \mathbf{m}_p
```
其中
```math
a_J = \frac{\hbar \gamma}{2 \mu_0 e} \frac{P J}{t M_s}, \qquad b_J = \xi a_J 
```

## 类阻尼与类场力矩

LLG 方程可以扩展以包含类阻尼（damping-like）与类场（field-like）力矩，这对模拟自旋轨道矩（SOT）或 $\epsilon(\theta)$ 为常数的 Slonczewski 力矩至关重要。修改后的 LLG 方程为：

```math
\frac{\partial \mathbf{m}}{\partial t} = -\gamma \, \mathbf{m} \times \mathbf{H}_{\text{eff}} 
+ \alpha \, \mathbf{m} \times \frac{\partial \mathbf{m}}{\partial t} 
- a_J \, \mathbf{m} \times (\mathbf{m} \times \mathbf{p}) 
- b_J \, \mathbf{m} \times \mathbf{p}
```

其中 $a_J$ 和 $b_J$ 分别是类阻尼与类场力矩的系数，$\mathbf{p}$ 是单位矢量。特别是在 SOT 的语境下，$\mathbf{p}$ 表示由自旋轨道耦合产生的自旋流方向。

## 有效场表述

在数值模拟中，自旋转移矩（STT）可以等价地表示为有效场贡献。把 STT 表示为有效场的 LLG 方程为：

```math
\frac{\partial \mathbf{m}}{\partial t} = -\gamma \, \mathbf{m} \times (\mathbf{H}_{\text{eff}} + \mathbf{H}_{\text{stt}}) 
+ \alpha \, \mathbf{m} \times \frac{\partial \mathbf{m}}{\partial t}
```

在 Zhang-Li 模型中，有效场为
```math
\mathbf{H}_{\text{stt}} = -\frac{1}{\gamma} \left[ \mathbf{m} \times (\mathbf{u} \cdot \nabla) \mathbf{m} + \beta \, (\mathbf{u} \cdot \nabla) \mathbf{m} \right]
```

对于类阻尼与类场力矩，有效场可以导出为：
```math
\mathbf{H}_{\text{stt}} = \frac{1}{\gamma} \left( a_J \, \mathbf{m} \times \mathbf{p} + b_J \, \mathbf{p} \right)
```

## SLLG 方程

SLLG 方程，即包含随机场 $\mathbf{b}$ 的 LLG 方程，为

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times (\mathbf{H}_\mathrm{eff} +\mathbf{b}) + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
```
热涨落假设为高斯白噪声，即热噪声 $\mathbf{b}$ 满足

```math
\left< \mathbf{b} \right> = 0, \;\;\; \left< \mathbf{b}_i^u \cdot \mathbf{b}_j^v \right> = 2 D \delta_{ij} \delta_{uv}
```

其中 $i$ 和 $j$ 是笛卡尔指标，$u$ 和 $v$ 标记磁化分量，$\left< \cdot , \cdot \right>$
表示对不同涨落场实现的平均。并且

```math
D = \frac{\alpha k_B T}{\gamma \mu_s}.
```

对于微磁情形，$D$ 为

```math
D = \frac{\alpha k_B T}{\mu_0 M_s \gamma \Delta V}.
```
这等价于随机场

```math
\mathbf{b}^u = \eta \sqrt \frac{2 \alpha k_B T}{\mu_0 M_s \gamma \Delta V dt}
```

其中 $\eta$ 是服从正态分布的随机数。

## 最速下降法

我们为复杂系统提供了最速下降能量极小化方法，其形式为

```math
x_{k+1} = x_k + \alpha_k d_k
```

其中

```math
d_k = - \nabla f(x_k)
```

对于微磁学，有

```math
\mathbf{m}_{k+1} = \mathbf{m}_{k} - {\tau}_k \mathbf{m}_k  \times (\mathbf{m}_k \times \mathbf{H}_{\mathrm{eff}})
```

实际中我们使用以下更新规则来保持磁化矢量的模长不变。

```math
\boldsymbol{m}_{k+1}=\boldsymbol{m}_{k}-{\tau}_k \frac{\boldsymbol{m}_{k}+\boldsymbol{m}_{k+1}}{2} \times\left(\boldsymbol{m}_{k} \times \boldsymbol{H}_{\mathrm{eff}}\left(\boldsymbol{m}_{k}\right)\right)
```

```math
\boldsymbol{m}_{k+1}^2 = \boldsymbol{m}_{k}^2
```

由该方程可得：

```math
(1+\frac{{\tau}_k^2}{4} \boldsymbol{f}_k^2)\mathbf{m}_{k+1} =
(1-\frac{{\tau}_k^2}{4} \boldsymbol{f}_k^2)\mathbf{m}_{k} -  {\tau}_k \mathbf{g}_k
```

其中

```math
\begin{aligned}
\mathbf{f}_k& = \mathbf{m}_k \times \mathbf{H}_{\mathrm{eff}}
\\\boldsymbol{g}_{k} &=\boldsymbol{m}_{k} \times\left(\boldsymbol{m}_{k} \times \boldsymbol{H}_{\mathrm{eff}}\right)
\end{aligned}
```

步长 $\tau_k$ 可以计算为

```math
\tau_{k}^{1}=\frac{\sum_{i} \boldsymbol{s}_{k-1}^{i} \cdot \boldsymbol{s}_{k-1}^{i}}{\sum_{i} \boldsymbol{s}_{k-1}^{i} \cdot \boldsymbol{y}_{k-1}^{i}} \quad, \quad \tau_{k}^{2}=\frac{\sum_{i} \boldsymbol{s}_{k-1}^{i} \cdot \boldsymbol{y}_{k-1}^{i}}{\sum_{i} \boldsymbol{y}_{k-1}^{i} \cdot \boldsymbol{y}_{k-1}^{i}}
```

其中

```math
\begin{aligned}  \boldsymbol{s}_{k-1} &=\boldsymbol{m}_{k}-\boldsymbol{m}_{k-1} \\ \boldsymbol{y}_{k-1} &=\boldsymbol{g}_{k}-\boldsymbol{g}_{k-1} \end{aligned}
```

## 蒙特卡洛模拟

实现的能量为

```math
\mathcal{H} = -\sum_{\langle i, j\rangle} \left( J_x S_i^x  S_j^x + J_y S_i^y  S_j^y + J_z S_i^z  S_j^z \right)

+ \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{S}_{i} \times \mathbf{S}_{j}\right)

- K \sum_{i}\left(\mathbf{u} \cdot \mathbf{S}_i\right)^2 - \sum_{i} \mathbf{H} \cdot \mathbf{S}_i
```

其中 $\mathbf{S}_i$ 是格点 _i_ 上经典自旋的单位矢量。

对于界面 DMI，

```math
\mathbf{D}_{i j} = D \hat{z} \times \hat{r}_{ij}  + D_z^{j} \hat{z}
```
而对于体 DMI，

```math
\mathbf{D}_{i j} = D \hat{r}_{ij}
```

注意蒙特卡洛只支持三角和立方网格。

## NEB（弹性带方法）

NEB 是一种寻找两个状态之间 MEP（最小能量路径）的链式方法。首先需要构建一条包含若干 image（每个 image 是磁化的一份拷贝）的链，然后对体系做弛豫。与初态和末态对应的两个端点 image 会被固定，因为它们是用户给定的能量态。由所有自由 image 组成的体系会被弛豫以降低总能量，这与禁用进动项时用 LLG 方程弛豫磁体系非常相似。一个重要的区别是：LLG 方程中的有效场是体系能量对磁化的泛函导数，而 NEB 中 image _n_ 的有效场还要包含其近邻（即 image _n-1_ 和 _n+1_）的影响。这种影响由所谓的切线（tangent）描述：弛豫整个体系时只保留有效场垂直于切线的分量。

假设整个体系有 _N_ 个 image

```math
\mathbf{Z} = [\mathbf{Y}_1, \mathbf{Y}_2, ..., \mathbf{Y}_N]
```

其中

```math
\mathbf{Y}_i = [m_{1x}, m_{1y}, m_{1z}, ..., m_{nx}, m_{ny}, m_{nz}]
```

每个 image 有 _n_ 个自旋。要弛豫体系，可以求解方程

```math
\frac{\mathbf{Y}_i}{\partial t} = - \mathbf{Y}_i \times (\mathbf{Y}_i \times \mathbf{G}_i)
```

其中 $\mathbf{G}_i$ 是有效场，可以计算为

```math
\mathbf{Y}_i = \mathbf{H}_i - (\mathbf{H}_i \cdot \mathbf{t}_i) \mathbf{t}_i +  \mathbf{F}_i
```

$\mathbf{H}_i$ 是常规的微磁有效场，$\mathbf{t}_i$ 是切线，$\mathbf{F}_i$ 是用于调节 image 间距的力。

```math
 \mathbf{F}_i = k (|\mathbf{Y}_{i+1}-\mathbf{Y}_{i}|-|\mathbf{Y}_{i}-\mathbf{Y}_{i-1}|) \mathbf{t}_i
```

image $\mathbf{Y}_{i}$ 与 $\mathbf{Y}_{j}$ 之间的距离定义为

```math
 L = \left [ \sum_k (L_k^{i,j})^2 \right ] ^{1/2}
```

其中 $L_k^{i,j}$ 是点 `k` 的测地距离，可以用 Vincenty 公式计算。

切线可以按如下方式计算

```math
\mathbf{t}_i^+ =  \mathbf{Y}_{i+1}-\mathbf{Y}_{i}\\
\mathbf{t}_i^- =  \mathbf{Y}_{i}-\mathbf{Y}_{i-1}
```

详细的方程见 [Journal of Chemical Physics 113, 22 (2000)] 和 [Computer Physics Communications 196 (2015) 335–347]。
