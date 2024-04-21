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

The implemented energy reads

```math
\mathcal{H} = -\sum_{\langle i, j\rangle} \left( J_x S_i^x  S_j^x + J_y S_i^y  S_j^y + J_z S_i^z  S_j^z \right)

+ \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{S}_{i} \times \mathbf{S}_{j}\right)

- K \sum_{i}\left(\mathbf{u} \cdot \mathbf{S}_i\right)^2 - \sum_{i} \mathbf{H} \cdot \mathbf{S}_i
```

where $\mathbf{S}_i$ is unit vector of the classical spin at site _i_.

For interfacial DMI,

```math
\mathbf{D}_{i j} = D \hat{z} \times \hat{r}_{ij}  + D_z^{j} \hat{z}
```
while for Bulk DMI, 

```math
\mathbf{D}_{i j} = D \hat{r}_{ij}
```

Note that the Monte Carlo only works for triangular and cubic meshes.

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
