# Equations

## Atomistic spin model

The basic assumption of the atomistic spin model is that each lattice site is associated with a magnetic moment $\mu_s$. For metal systems with quenched orbital moments, the magnetic moment is mainly related to its spin angular momentum
```math
\mathbf{\mu} = - g \mu_B \mathbf{S} = - \hbar  \gamma \mathbf{S}
```
where $\mu_B=e \hbar /(2m)$ is the Bohr magneton, $e(>0)$ is the electron charge, 
$\gamma=g\mu_B/\hbar (>0) $ is the gyromagnetic ratio, $g=2$ is the g-factor. 
The LLG equation governs the dynamics of the magnetic moment, 
which reads

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H}_\mathrm{eff} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
```

where $\mathbf{m}$ is the unit vector of the magnetic moment, $\mathbf{H}_\mathrm{eff}$ is the effective field.

Unlike the continuous micromagnetic model, the atomistic spin model is discrete, and thus, the effective field $\mathbf{H}_\mathrm{eff}$ is defined as the partial derivative of the total Hamiltonian with respect to $\mathbf{m}$, i.e., 

```math
\mathbf{H}_{\mathrm{eff}}=-\frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \mathbf{m}}
```

where $\mathcal{H}$ is the total Hamiltonian including the exchange interaction, Dzyaloshinskii-Moriya interaction, dipolar interaction, anisotropy interaction,  Zeeman interaction and so on.  


The exchange interaction is given by

```math
\mathcal{H}_\mathrm{ex} = -J \sum_{\langle i, j\rangle} \mathbf{m}_{i} \cdot \mathbf{m}_{j}
```
where $J$ denotes the exchange constant and $\langle i,j \rangle$ represents a unique pair between lattice sites $i$ and $j$ and we assume that the summation is taken only once for each pair. A positive $J$ results in the ferromagnetic state while a negative $J$ leads to the antiferromagnetic state, which can be seen by minimizing the exchange Hamiltonian. To solve the LLG equation numerically, we need the effective field at each site, which is given 

```math
\mathbf{H}_\mathrm{ex, i} = \frac{J}{\mu_s} \sum_{\langle i, j\rangle} \mathbf{m}_{j}.
```

Similarly, the Dzyaloshinskii-Moriya interaction reads

```math
\mathcal{H}_\mathrm{dmi} =  \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right)
```
where $\mathbf{D}_{i j}$ is the DM vector.  Typically, there are two situations: (1) $\mathbf{D}_{i j} = D \hat{r}_{ij}$, which corresponds to the Bulk DMI. (2) $\mathbf{D}_{i j} = D \hat{r}_{ij} \times \hat{z}$, which corresponds to the interfacial DMI. Unlike the exchange interaction, the DM Hamiltonian is minimal when the two adjacent magnetic moments are perpendicular. The corresponding effective field can be computed by
```math
\mathbf{H}_\mathrm{dmi, i} = \frac{1}{\mu_s} \sum_{\langle i, j\rangle} \mathbf{D}_{i j} \times \mathbf{m}_{j}.
```

The other two common interactions are the anisotropy and the Zeeman, which are given by

```math
\mathcal{H}_\mathrm{an} = - K \sum_{i}\left(\mathbf{m}_{i} \cdot \hat{u}\right)^{2} \\
\mathcal{H}_\mathrm{ze} =  - \mu_s \mathbf{m}_i \cdot \mathbf{H}
```

The effective field of the anisotropy is 
```math
\mathbf{H}_\mathrm{an, i} = \frac{2K}{\mu_s}  (\mathbf{m}_{i} \cdot \hat{u}) \hat{u}.
```

### Multiferroic Insulators

For multiferroics such as Cu$_2$OSeO$_3$, the noncollinear spin texture induces a local electric polarization via the the $d-p$ hybridization mechanism. The local electric dipole moment $\mathbf{P}$ depends on the direction of the applied static magnetic field relative to the crystallographic axes.[PRL 128, 037201 (2022)] 


##### $H_0//001$:
```math
\mathbf{P}_i=\lambda\left(-m_{i, z} m_{i, x}, m_{i, y} m_{i, z}, \frac{-m_{i, x}^2+m_{i, y}^2}{2}\right)
```

##### $H_0//110$:
```math
\mathbf{P}_i=\lambda\left(-m_{i, x} m_{i, y}, \frac{-m_{i, x}^2+m_{i, z}^2}{2}, m_{i, y} m_{i, z}\right)
```

##### $H_0//111$:

```math
\mathbf{P}_{i,x}=-\lambda \frac{m_{i, x}\left(\sqrt{2} m_{i, y}+m_{i, z}\right)}{\sqrt{3}} \\

\mathbf{P}_{i,y} = \lambda \frac{-m_{i, x}^2+m_{i, y}\left(m_{i, y}-\sqrt{2} m_{i, z}\right)}{\sqrt{6}} \\

\mathbf{P}_{i,z} =
-\lambda \frac{m_{i, x}^2+m_{i, y}^2-2 m_{i, z}^2}{2 \sqrt{3}}
```

The Hamiltonian related to the high-frequency lasers is given by 

```math
\mathcal{H}_\mathrm{laser} =  -\sum_{i} \mu_s \mathbf{m}_i \cdot \mathbf{H}(t) - \sum_{i} \mathbf{P}_i \cdot \mathbf{E}(t)
```

The corresponding effective fields associated with the electric field are

##### $H_0//001$:
```math
H_x = \frac{\lambda}{\mu_s}  (-E_z m_x  - E_x m_z )\\
H_y = \frac{\lambda}{\mu_s}  (E_z m_y  + E_y m_z )\\
H_z = \frac{\lambda}{\mu_s}  (-E_x m_x  + E_y m_y )\\
```

##### $H_0//110$:
```math
H_x = \frac{\lambda}{\mu_s}  (-E_y m_x  - E_x m_y )\\
H_y = \frac{\lambda}{\mu_s}  (-E_x m_x  + E_z m_z )\\
H_z = \frac{\lambda}{\mu_s}  (E_z m_y  + E_y m_z )\\
```

##### $H_0//111$:
```math
H_x = -\frac{\lambda}{\sqrt{3}\mu_s}  (\sqrt{2}(E_y m_x + E_x m_y)  +E_z m_x + E_x m_z )\\
H_y = -\frac{\lambda}{\sqrt{3}\mu_s}  (\sqrt{2}(E_x m_x - E_y m_y)  +E_z m_y + E_y m_z )\\
H_z = -\frac{\lambda}{\sqrt{3}\mu_s}  (E_x m_x  + E_y m_y - 2 E_z m_z)\\
```

## Micromagnetic model

In micromagnetics, the effective field can be computed from the total micromagnetic energy

```math
\mathbf{H}_{\mathrm{eff}}=-\frac{1}{\mu_{0} M_{s}} \frac{\delta E}{\delta \mathbf{m}}
```

The typical energy terms are

- **Exchange energy**

```math
  E_\mathrm{ex} = \int_{V} A (\nabla \mathbf{m})^2 \mathrm{d}V
```

  where $(\nabla \mathbf{m})^{2}=\left(\nabla m_{x}\right)^{2}+\left(\nabla m_{y}\right)^{2}+\left(\nabla m_{z}\right)^{2}$. So the corresponding effective field is

```math
  \mathbf{H}_{\mathrm{ex}}=\frac{2 A}{\mu_{0} M_{s}} \nabla^{2} \mathbf{m}
```

- **Zeeman energy**

```math
  E_\mathrm{ex} = -  \mu_0 \int_{V}  \mathbf{H} \cdot \mathbf{M} \mathrm{d}V
```

  as expected, the effective field is $\mathbf{H}$.

- **Anisotropy**

  The uniaxial anisotropy energy is given by

```math
  E_\mathrm{anis} = -\int_{V} K_{u} (\mathbf{m} \cdot \hat{u})^2 \, dV
```

  from which the effective field can be computed as

```math
  \mathbf{H}_{\mathrm{an}}=\frac{2 K_u}{\mu_0 M_s}\left(\mathbf{m} \cdot \hat{u}\right) \hat{u}
```

- **Cubic Anisotropy**

  The cubic anisotropy energy is given by

```math
  E_\mathrm{cubic} = -\int_{V} K_c (m_x^4 + m_y^4 + m_z^4) \, dV
```

  and thus the corresponding effective field reads

```math
  \mathbf{H}_{\mathrm{cubic}}= \frac{4 K_c}{\mu_0 M_s}  
  \left( m_x^3 \mathbf{e}_x + m_y^3 \mathbf{e}_y + m_z^3 \mathbf{e}_z \right)
```

- **Bulk DMI energy** The Bulk DMI energy reads

```math
  E_{\mathrm{dmi}} = \int_V D \mathbf{m} \cdot (\nabla \times \mathbf{m}) \, \mathrm{d}V
```

  so the effective field is

```math
  \mathbf{H}_\mathrm{D}=-\frac{2 D}{\mu_{0} M_{s}}(\nabla \times \mathbf{m})
```

- **Magnetostatic energy**

```math
  E_{\mathrm{d}}=-\frac{\mu_{0}}{2} \int_{V} \mathbf{H}_{\mathrm{d}}(\mathbf{r}) \cdot \mathbf{M}(\mathbf{r}) d V
```

```math
  \mathbf{H}_{\mathrm{d}}(\mathbf{r})=\frac{1}{4 \pi}\left(\int_{V} \rho_{m}\left(\mathbf{r}^{\prime}\right) \frac{\mathbf{r}-\mathbf{r}^{\prime}}{\left|\mathbf{r}-\mathbf{r}^{\prime}\right|^{3}} \mathrm{d}^{3} r^{\prime}+\int_{S} \sigma_{m}\left(\mathbf{r}^{\prime}\right) \frac{\mathbf{r}-\mathbf{r}^{\prime}}{\left|\mathbf{r}-\mathbf{r}^{\prime}\right|^{3}} \mathrm{d}^{2} r^{\prime}\right)
```

## LLG equation

The driver `LLG` solves the standard LLG equation, which can be written as

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
```

and the corresponding LL form is given by

```math
(1+\alpha^2)\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} - \alpha \gamma \mathbf{m} \times (\mathbf{m} \times \mathbf{H})
```

## LLG equation with extensions

### Spin transfer torque

In micromagnetics, the spin transfer torque is modelled with the extended LLG equation Zhang-Li extension, 
which reads

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
+ (\mathbf{u} \cdot \nabla) \mathbf{m} - \beta [\mathbf{m}\times (\mathbf{u} \cdot \nabla)\mathbf{m}]
```
where 
```math
\mathbf{u} = \frac{p g \mu_B}{2 e M_s} \mathbf{j}
```
represents the strength of the current. The unit of $\mathbf{u}$ is m/s. 
In the definition of  $\mathbf{u}$,
$p$ is the spin polarization of the electric current, $e(>0)$ is the elementary charge,
$M_s$ is the saturation magnetization and $\mathbf{j}$ is electric current density.

For the atomistic model, 
```math
\mathbf{u} = \frac{p g \mu_B a^3}{2 e \mu_s} \mathbf{j} = \frac{p a^3}{2 e S} \mathbf{j}.
```
where $S = |\mathbf{S}|$ is the length of local spin.

In JuMag, the the extended LLG equation Zhang-Li extension is implemented in the driver `LLG_STT`.  
Moreover, the driver `LLG_CPP` implements the LLG equation with spin transfer torque
for the current-perpendicular-to-plane (CPP) case,

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t} - a_J \mathbf{m} \times (\mathbf{m} \times \mathbf{p})
 -  b_J \mathbf{m} \times \mathbf{p}
```
The spin valve structures and  spin orbit torques can use the `LLG_CPP` driver.

The driver `LLG_STT_CPP` has put them together,

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
+ (\mathbf{u} \cdot \nabla) \mathbf{m} - \beta [\mathbf{m}\times (\mathbf{u} \cdot \nabla)\mathbf{m}] - a_J \mathbf{m} \times (\mathbf{m} \times \mathbf{p})
 -  b_J \mathbf{m} \times \mathbf{p}
```

## SLLG equation

The SLLG equation, i.e., LLG equation including the stochastic field $\mathbf{b}$, is given by

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times (\mathbf{H}_\mathrm{eff} +\mathbf{b}) + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
```
The thermal fluctuation is assumed to be a Gaussian white noise, i.e., the thermal noise $\mathbf{b}$ obeys the  properties

```math
\left< \mathbf{b} \right> = 0, \;\;\; \left< \mathbf{b}_i^u \cdot \mathbf{b}_j^v \right> = 2 D \delta_{ij} \delta_{uv}
```

where $i$ and $j$ are Cartesian indices, $u$ and $v$ indicate the magnetization components and $\left< \cdot , \cdot \right>$
represents the average taken over different realizations of the fluctuating field. And

```math
D = \frac{\alpha k_B T}{\gamma \mu_s}.
```

For the micromagnetic case, $D$ is given as

```math
D = \frac{\alpha k_B T}{\mu_0 M_s \gamma \Delta V}.
```
which is equivalent to a stochastic field  

```math
\mathbf{b}^u = \eta \sqrt \frac{2 \alpha k_B T}{\mu_0 M_s \gamma \Delta V dt}
```

where $\eta$ is a random number follows the normal distribution.

## Steepest descent method

We provide a steepest descent energy minimization method for a complicated system, which is of the form

```math
x_{k+1} = x_k + \alpha_k d_k
```

where

```math
d_k = - \nabla f(x_k)
```

And for the micromagnetics, we have

```math
\mathbf{m}_{k+1} = \mathbf{m}_{k} - {\tau}_k \mathbf{m}_k  \times (\mathbf{m}_k \times \mathbf{H}_{\mathrm{eff}})
```

In practice, we use the following update rule to keep the magnetization vector normalized.

```math
\boldsymbol{m}_{k+1}=\boldsymbol{m}_{k}-{\tau}_k \frac{\boldsymbol{m}_{k}+\boldsymbol{m}_{k+1}}{2} \times\left(\boldsymbol{m}_{k} \times \boldsymbol{H}_{\mathrm{eff}}\left(\boldsymbol{m}_{k}\right)\right)
```

```math
\boldsymbol{m}_{k+1}^2 = \boldsymbol{m}_{k}^2
```

From the equation we have:

```math
(1+\frac{{\tau}_k^2}{4} \boldsymbol{f}_k^2)\mathbf{m}_{k+1} =
(1-\frac{{\tau}_k^2}{4} \boldsymbol{f}_k^2)\mathbf{m}_{k} -  {\tau}_k \mathbf{g}_k
```

where

```math
\begin{aligned}
\mathbf{f}_k& = \mathbf{m}_k \times \mathbf{H}_{\mathrm{eff}}
\\\boldsymbol{g}_{k} &=\boldsymbol{m}_{k} \times\left(\boldsymbol{m}_{k} \times \boldsymbol{H}_{\mathrm{eff}}\right)
\end{aligned}
```

The step size $\tau_k$ can be computed by

```math
\tau_{k}^{1}=\frac{\sum_{i} \boldsymbol{s}_{k-1}^{i} \cdot \boldsymbol{s}_{k-1}^{i}}{\sum_{i} \boldsymbol{s}_{k-1}^{i} \cdot \boldsymbol{y}_{k-1}^{i}} \quad, \quad \tau_{k}^{2}=\frac{\sum_{i} \boldsymbol{s}_{k-1}^{i} \cdot \boldsymbol{y}_{k-1}^{i}}{\sum_{i} \boldsymbol{y}_{k-1}^{i} \cdot \boldsymbol{y}_{k-1}^{i}}
```

where

```math
\begin{aligned}  \boldsymbol{s}_{k-1} &=\boldsymbol{m}_{k}-\boldsymbol{m}_{k-1} \\ \boldsymbol{y}_{k-1} &=\boldsymbol{g}_{k}-\boldsymbol{g}_{k-1} \end{aligned}
```

## Monte Carlo Simulation

For cubic mesh, the implemented energy reads

```math
\mathcal{H} = - J \sum_{\langle \langle i, j\rangle \rangle} \mathbf{m}_{i} \cdot \mathbf{m}_{j}
+\sum_{\langle \langle i, j \rangle \rangle }  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right)
- K \sum_{i}\left(\mathbf{u} \cdot \mathbf{m}_i\right)^2 - \sum_{i} \mathbf{H} \cdot \mathbf{m}_i
```

where $\mathbf{m}_i$ is unit vector of the classical spin at site _i_.

For triangular mesh (2D), the system energy reads

```math
H= \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{S}_{i} \times \mathbf{S}_{j}\right)
-J \sum_{\langle i, j\rangle} \mathbf{S}_{i} \cdot \mathbf{S}_{j}- \lambda \sum_{\langle i, j\rangle} S_{i}^{z} S_{j}^{z}
-K \sum_{i}\left(S_{i}^{z}\right)^{2}
```

where

```math
\mathbf{D}_{i j} = D \hat{z} \times \hat{r}_{ij}  + D_z^{j} \hat{z}
```

## NEB (Nudged elastic band)

NEB is a chain method to find the MEP (minimum energy path) between two states. To start, we need to construct a chain including several images (each image is a copy of the magnetization) and then relax the system. Two ends images that corresponding to the initial and final states will be pinned as they are the energy states given by the users. The system contain all free images will be relaxed to reduce the total energy, which is very similar to the case that relaxing the magnetic system using LLG equation if one disables the precession term. One significant difference is that the effective field in LLG equation is the functional derivative of the system energy with respect to magnetization while in NEB the effective field of image _n_ should also contain the influence of its neighbours (i.e., the images _n-1_ and _n+1_). This influence is described by the so-called tangents: only the perpendicaular part of the effective field is kept when relaxing the whole system.

Assume that the whole system has _N_ images

```math
\mathbf{Z} = [\mathbf{Y}_1, \mathbf{Y}_2, ..., \mathbf{Y}_N]
```

where

```math
\mathbf{Y}_i = [m_{1x}, m_{1y}, m_{1z}, ..., m_{nx}, m_{ny}, m_{nz}]
```

each image has _n_ spins. To relax the system, we could solve the equation

```math
\frac{\mathbf{Y}_i}{\partial t} = - \mathbf{Y}_i \times (\mathbf{Y}_i \times \mathbf{G}_i)
```

where $\mathbf{G}_i$ is effective field that can be computed as

```math
\mathbf{Y}_i = \mathbf{H}_i - (\mathbf{H}_i \cdot \mathbf{t}_i) \mathbf{t}_i +  \mathbf{F}_i
```

The $\mathbf{H}_i$ is the normal micromagnetic effective field, $\mathbf{t}_i$ is the tangent and $\mathbf{F}_i$ is a force that can be used to adjust the distance between images.

```math
 \mathbf{F}_i = k (|\mathbf{Y}_{i+1}-\mathbf{Y}_{i}|-|\mathbf{Y}_{i}-\mathbf{Y}_{i-1}|) \mathbf{t}_i
```

The distance bewteen images $\mathbf{Y}_{i}$ and $\mathbf{Y}_{j}$ is defined as

```math
 L = \left [ \sum_k (L_k^{i,j})^2 \right ] ^{1/2}
```

where $L_k^{i,j}$ is the geodesic distance of point `k` that can be computed using Vincenty's formula.

The tangents can be computed as follows

```math
\mathbf{t}_i^+ =  \mathbf{Y}_{i+1}-\mathbf{Y}_{i}\\
\mathbf{t}_i^- =  \mathbf{Y}_{i}-\mathbf{Y}_{i-1}
```

The detailed equations can be found @ [Journal of Chemical Physics 113, 22 (2000)] and [Computer Physics Communications 196 (2015) 335–347].

## Eigenvalue Method

Micromagnetic system and the corresponding atomistic model are the classical system. The resonance frequencies and the spatial resonance modes can be obtained using the eigenvalue method. In JuMag, we implemented a simple eigenvalue method for the following Hamiltonian of the atomistic model

```math
\mathcal{H} = -J \sum_{\langle i, j\rangle} \mathbf{m}_{i} \cdot \mathbf{m}_{j} + \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right) - K \sum_{i}\left(m_{i}^{z}\right)^{2} - \mu_s \mathbf{m}_i \cdot \mathbf{H}
```

For a given ground state or metastable spin configuration

```math
\mathbf{m}_0=(\sin\theta \cos \phi, \sin\theta \sin \phi, \cos\theta)^T,
```

one can construct a local coordinate system such that

```math
\mathbf{m} = w \mathbf{m}_0 + u \mathbf{e}_\theta + v  \mathbf{e}_\phi,
```

where $\mathbf{e}_\theta =(\cos \theta \cos \phi, \cos\theta \sin\phi,-\sin \theta)^T$ and $\mathbf{e}_\phi=(-\sin \phi, \cos\phi,0)^T$. In matrix form,

```math
\begin{bmatrix}
m_x \\ m_y \\ m_z
 \end{bmatrix} =
 \begin{bmatrix}
 \cos \phi \cos \theta & -\sin \phi & \sin \theta \cos \phi \\
 \sin \phi \cos \theta &  \cos \phi & \sin \theta \sin \phi \\
 -\sin \theta & 0 &  \cos \theta \\
 \end{bmatrix}
  \begin{bmatrix}
 u \\ v \\ w
 \end{bmatrix}.
```

Under this transformation, the unperturbed spin configuration $\mathbf{m}_0$ corresponds to $u=0, v=0$ and $w=1$. In this local coordinate, the effective fields are given by

```math
\mathbf{H}_\mathrm{eff} = - \frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \mathbf{m}}
= - \frac{1}{\mu_s}\left(  \frac{\partial \mathcal{H}}{\partial w} \mathbf{m}_0 +
\frac{\partial \mathcal{H}}{\partial u} \mathbf{e}_\theta + \frac{\partial \mathcal{H}}{\partial v} \mathbf{e}_\phi \right).
```

Substituting the effective fields into the LLG equation, we obtain for the case that $\alpha=0$

```math
\begin{matrix}
\dot{u} = -\gamma (- w H_v + v H_w),\\
\dot{v} = -\gamma ( w H_u - u H_w),
\end{matrix}
```

where $H_w$, $H_u$ and $H_v$ are effective fields along $\mathbf{m}_0$, $\mathbf{e}_\theta$ and $\mathbf{e}_\phi$, respectively. To linearize the LLG equation, we assume $|u| \ll 1$, $|v| \ll 1$ and thus $w \approx 1- (1/2)(u^2+v^2)$. Moreover, we look for the solutions such that $u=\tilde{u}e^{-i\omega t}$ and $v=\tilde{v}e^{-i\omega t}$. Therefore, we arrive at

```math
\begin{matrix}
 -i \omega \mu_s  \tilde{u}  =  \gamma (\tilde{h}_v - \tilde{v} \tilde{H}_w) , \\
-i \omega \mu_s \tilde{v}  =  - \gamma (\tilde{h}_u-\tilde{u} \tilde{H}_w),
\end{matrix}
```

where we have ignored the higher-order terms. In addition, the above equation can be rewritten in a matrix form~\cite{Lin2014a}

```math
\frac{i \omega \mu_s}{\gamma}  
\begin{bmatrix}
\mathbf{u} \\ \mathbf{v}
 \end{bmatrix} = \mathbf{A}  
 \begin{bmatrix}
\mathbf{u} \\ \mathbf{v}
 \end{bmatrix},
```

where we have introduced two vectors $\mathbf{u}=(\tilde{u}_1, \tilde{u}_2, ..., \tilde{u}_n)^T$ and $\mathbf{v}=(\tilde{v}_1, \tilde{v}_2, ..., \tilde{v}_n)^T$. Therefore, the normal modes of the system can be obtained through solving the eigenvalues of the matrix $\mathbf{A}$. The eigenvalues are pure imaginary numbers since $\omega$ is real. The related effective fields are given by

```math
\begin{aligned}
\tilde{H}_w = -(\partial \mathcal{H}^{(0)}/{\partial w}) \big|_{w=1} \\
\tilde{h}_u = -(\partial \mathcal{H}^{(2)}/{\partial u}) \big|_{w=1, u=\tilde{u}, v=\tilde{v}} \\
\tilde{h}_v = -(\partial \mathcal{H}^{(2)}/{\partial v}) \big|_{w=1, u=\tilde{u}, v=\tilde{v}} .
\end{aligned}
```

- **Exchange interaction**

```math
  \begin{aligned}
  \tilde{H}^\mathrm{ex}_{w, i} = J  \sum_j \left [ \cos \theta_i \cos \theta_j + \sin \theta_i \sin \theta_j  \cos(\phi_i-\phi_j) \right],\\
  \tilde{h}^\mathrm{ex}_{u, i} =  J \sum_j \left [ \tilde{u}_j \cos(\phi_i-\phi_j) \cos \theta_i \cos \theta_j + \tilde{u}_j \sin \theta_i \sin \theta_j + \tilde{v}_j \cos\theta_i \sin (\phi_i - \phi_j)  \right], \\
  \tilde{h}^\mathrm{an}_{v, i} =   J \sum_j \left [  -\tilde{u}_j \cos\theta_j \sin (\phi_i - \phi_j) + \tilde{v}_j \cos(\phi_i-\phi_j)  \right].
  \end{aligned}
```

- **DMI**

```math
  \begin{aligned}
  \tilde{H}^\mathrm{dmi}_{w, i} =  \sum_{j \in X} D_{ij}  \left(  \sin \theta_j \sin \phi_j \cos \theta_i - \sin \theta_i \sin \phi_i \cos \theta_j \right) + \\
     \sum_{j \in Y} D_{ij}  \left(  \sin \theta_i \cos \phi_i \cos \theta_j - \sin \theta_j \cos \phi_j \cos \theta_i  \right) + \\
     \sum_{j \in Z} D_{ij}   \left [  \sin \theta_i \sin \theta_j \sin(\phi_i-\phi_j) \right],
  \end{aligned}
```

  where $D_{ij}=D\mathrm{sgn}(j-i)$ with $\mathrm{sgn}(x)$ the sign function. The sets $X$, $Y$ and $Z$ represent the neighbours of site $i$ in $x$-, $y$- and $z$-axis. Meanwhile, $\tilde{h}^\mathrm{dmi}_{u, i}$ and $\tilde{h}^\mathrm{dmi}_{v, i}$ are given by

```math
  \begin{aligned}
  \tilde{h}^\mathrm{dmi}_{u, i} =  \sum_{j \in X} D_{ij} \left(
  \tilde{u}_j \sin \theta_j \sin \phi_i \cos \theta_i - \tilde{u}_j \sin \theta_i \sin \phi_j \cos \theta_j - \tilde{v}_j \sin \theta_i \cos\phi_j \right) +\\
  \sum_{j\in Y} D_{ij} \left (\tilde{u}_j \sin \theta_i \cos \phi_j \cos \theta_j - \tilde{u}_j \sin \theta_j \cos \phi_i \cos \theta_i - \tilde{v}_j \sin \theta_i \sin\phi_j  \right) + \\
  \sum_{j\in Z} D_{ij} \cos \theta_i  \left [ \tilde{u}_j  \cos \theta_j \sin(\phi_i- \phi_j) -\tilde{v}_j \cos(\phi_i-\phi_j)  \right ],\\
  \tilde{h}^\mathrm{dmi}_{v, i} = \sum_{j \in X} D_{ij}  \tilde{u}_j \sin \theta_j \cos\phi_i + \textstyle \sum_{j \in Y} D_{ij} \tilde{v}_j \sin \theta_j  \sin\phi_i + \\
                  {\sum_{j \in Z}} D_{ij} \left[ \tilde{u}_j  \cos \theta_j  \cos(\phi_i -\phi_j) + \tilde{v}_j \sin (\phi_i-\phi_j) \right].
  \end{aligned}
```

- **Anisotropy** For anisotropies with $\mathcal{H}_{an} = - \sum_{i} (K_x m_{x,i}^2 + K_z m_{z,i}^2)$, these fields are given by

```math
  \begin{aligned}
  \tilde{H}^\mathrm{an}_{w, i} = 2K_x\cos^2 \phi_i \sin^2 \theta_i +  2K_z \cos^2 \theta_i, \\
  \tilde{h}^\mathrm{an}_{u, i} =  2K_x \cos \phi_i \cos \theta_i (\tilde{u}_i \cos\phi_i \cos \theta_i - \tilde{v}_i \sin\phi_i) +2 K_z \tilde{u}_i \sin^2\theta_i, \\
  \tilde{h}^\mathrm{an}_{v, i} = 2K_x \sin \phi_i (-\tilde{u}_i \cos\phi_i \cos \theta_i + \tilde{v}_i \sin\phi_i ).
  \end{aligned}
```

- **Zeeman Field** For a static external field $\mathbf{H}=(H_x, H_y, H_z)$, one obtains $\tilde{h}_u = \tilde{h}_v = 0$ and $\tilde{H}^a_{w, i} = H_z \cos \theta_i + H_x \cos \phi_i \sin \theta_i + H_y \sin \phi_i \sin \theta_i$.
