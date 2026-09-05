```@meta
ShareDefaultModule = true
```

# 典型磁结构的 LTEM 成像

本例对三种典型磁结构——**磁涡旋**与 **Bloch/Néel 两种螺旋度的 skyrmion**——
模拟洛伦兹透射电镜(LTEM)图像,并复现实验与 LTEM 模拟研究中的经典衬度特征:

* 涡旋的**圆形磁台阶相位**及其在散焦下的亮/暗环形衬度;
* skyrmion 衬度的**螺旋度依赖**:Bloch 型零倾角即可见,而纯 Néel 型零倾角
  *完全不可见*,只有倾斜样品后才出现;
* 欠焦与过焦之间的**衬度反转**。

导入模块并定义两个绘图辅助函数:

````@example ltem
using MicroMagnetic
using CairoMakie

function show_phase!(pos, phi, title)
    ax = Axis(pos; title=title, aspect=DataAspect())
    heatmap!(ax, phi; colormap=:balance)
    hidedecorations!(ax)
end

function show_image!(pos, img, title)
    ax = Axis(pos; title=title, aspect=DataAspect())
    heatmap!(ax, img; colormap=:grays, colorrange=(0.85, 1.15))
    hidedecorations!(ax)
end
````

本页使用 300 kV 显微镜,并令 `V0 = 0` 以跳过电相位(平均内电势),因此所有
衬度都是纯磁性的。散焦量单位 µm,倾角单位 rad。

## 圆盘中的磁涡旋

圆盘中的涡旋:面内磁化绕中心卷曲,中心为面外磁化核心。我们在 320 nm × 320 nm
窗口(胞元 2.5 nm)内构造半径 95 nm、厚度 24 nm 圆盘的标准涡旋构型:

````@example ltem
n = 128; dx = dy = 2.5e-9; dz = 24e-9
R = 95e-9
m = zeros(3, n, n, 1)
for i in 1:n, j in 1:n
    x = (i - 0.5) * dx - 64e-9
    y = (j - 0.5) * dy - 64e-9
    r = hypot(x, y)
    if r < R
        f = tanh(r / 20e-9)                      # 面内卷曲幅值
        m[1, i, j, 1] = -y / r * f
        m[2, i, j, 1] =  x / r * f
        m[3, i, j, 1] = -tanh(r / 8e-9)          # 向下的核心
    end
end
nothing #hide
````

涡旋的束流积分面内感应是方位角的,因此磁相位是一个**圆形台阶**:相位在核心
两侧各自恒定,在核心区域下降(本例厚度与磁化下约 −2.5 rad)。计算相位和三幅
Fresnel 散焦像:

````@example ltem
phi_v = compute_magnetic_phase(m; Ms=8e5, dx=dx, dy=dy, dz=dz)
im_v1 = LTEM(m; Ms=8e5, V=300, df=-300, V0=0, dx=dx, dy=dy, dz=dz)[2]
im_v2 = LTEM(m; Ms=8e5, V=300, df=+300, V0=0, dx=dx, dy=dy, dz=dz)[2]
im_v3 = LTEM(m; Ms=8e5, V=300, df=-300, V0=0, ty=0.30, dx=dx, dy=dy, dz=dz)[2]
nothing #hide
````

````@example ltem
fig = Figure(size=(1050, 760))
show_phase!(fig[1, 1], phi_v, "magnetic phase")
show_image!(fig[1, 2], im_v1, "df = -300 um, tilt 0")
show_image!(fig[2, 1], im_v2, "df = +300 um, tilt 0")
show_image!(fig[2, 2], im_v3, "df = -300 um, tilt 17 deg")
fig
````

要点:

* **相位图**显示以涡旋核心为中心的圆形磁台阶;
* 零倾角时散焦像呈现围绕核心位置的**同心亮/暗环**(台阶的 Fresnel 像),
  且欠焦与过焦的衬度**相反**;
* 倾斜样品会把面外核心混入束流积分感应,额外产生熟悉的**非对称核心衬度**。

## Skyrmion:Bloch 与 Néel 螺旋度

skyrmion 的面内感应绕核心以固定螺旋角 ``\chi`` 卷绕:
``\mathbf{m}_\perp \propto \cos\chi\,\hat{r} + \sin\chi\,\hat{\phi}``。
LTEM 相位核对投影感应的无旋(径向)部分完全不敏感,因此 **skyrmion 的零倾角
衬度按 ``\sin\chi`` 缩放**:Bloch skyrmion(``\chi = 90°``)零倾角即呈现教科书式
衬度;而纯 Néel skyrmion(``\chi = 0°``)零倾角相位*恒为零*,只有倾斜样品后才可见。

使用内置的 [`skyrmion`](@ref) 纹理(半径 30 nm,壁宽 15 nm),窗口 240 nm、
胞元 2 nm,膜厚 50 nm:

````@example ltem
n2 = 120; dx2 = dy2 = 2e-9; dz2 = 50e-9
function skyrmion_m(helicity)
    m = zeros(3, n2, n2, 1)
    sk = skyrmion(; center=(120e-9, 120e-9), R=30e-9, p=-1, v=1, phi=helicity)
    for i in 1:n2, j in 1:n2
        m[:, i, j, 1] .= sk((i - 0.5) * dx2, (j - 0.5) * dy2, 0.0)
    end
    return m
end
m_bloch = skyrmion_m(pi / 2)   # Bloch:面内分量为方位角的
m_neel  = skyrmion_m(0.0)      # Néel:面内分量为径向的
nothing #hide
````

````@example ltem
phi_b = compute_magnetic_phase(m_bloch; Ms=8e5, dx=dx2, dy=dy2, dz=dz2)
b_m = LTEM(m_bloch; Ms=8e5, V=300, df=-300, V0=0, dx=dx2, dy=dy2, dz=dz2)[2]
b_p = LTEM(m_bloch; Ms=8e5, V=300, df=+300, V0=0, dx=dx2, dy=dy2, dz=dz2)[2]
n_0 = LTEM(m_neel;  Ms=8e5, V=300, df=-300, V0=0, dx=dx2, dy=dy2, dz=dz2)[2]
n_t = LTEM(m_neel;  Ms=8e5, V=300, df=-300, V0=0, ty=0.30, dx=dx2, dy=dy2, dz=dz2)[2]
nothing #hide
````

````@example ltem
fig = Figure(size=(1250, 760))
show_phase!(fig[1, 1], phi_b, "Bloch: phase (tilt 0)")
show_image!(fig[1, 2], b_m, "Bloch: df = -300 um, tilt 0")
show_image!(fig[1, 3], b_p, "Bloch: df = +300 um, tilt 0")
show_phase!(fig[2, 1], phi_b .* 0, "Neel: phase (tilt 0) = 0")
show_image!(fig[2, 2], n_0, "Neel: df = -300 um, tilt 0")
show_image!(fig[2, 3], n_t, "Neel: df = -300 um, tilt 17 deg")
fig
````

读图:

* **Bloch skyrmion(上排)**:相位呈两瓣偶极结构,Fresnel 像是一对亮/暗瓣,
  欠焦与过焦时**亮暗互换**;
* **Néel skyrmion(下排)**:零倾角的相位与图像*完全平坦*——不倾角什么都
  看不到——而 17° 倾角下 skyrmion 周围出现清晰的环形衬度。

这一螺旋度依赖是洛伦兹显微术的著名结论,本例未引入任何可调参数即复现。

## 说明

* 零倾角的衬度完全来自面内感应;垂直磁化分量只有倾斜(或在样品边缘通过其
  1/r 长程尾部)才可见。
* 所有图像均以入射束强度归一化:亮表示电子多于背景,暗表示少于背景。
* 真实样品的电相位(`V0`)会在材料覆盖区上叠加一个强而无结构的背景。

## 参考文献

* J. N. Chapman and M. R. Scheinfein, "Transmission electron microscopies
  of magnetic microstructures", J. Magn. Magn. Mater. **200**, 729 (1999).
* C. Phatak, A. K. Petford-Long and M. De Graef, "Recent advances in
  Lorentz microscopy", Curr. Opin. Solid State Mater. Sci. **20**, 107
  (2016).
* S. K. Walton *et al.*, "MALTS: A tool to simulate Lorentz Transmission
  Electron Microscopy from micromagnetic simulations", J. Phys. D **46**,
  175005 (2013).
