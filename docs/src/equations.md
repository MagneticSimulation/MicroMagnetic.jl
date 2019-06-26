# Implemented equations

## Energies

- Exchange energy

  ```math
  E_\mathrm{ex} = \int_{V} A (\nabla \vec{m})^2 \mathrm{d}V
  ```

- Bulk DMI energy

  ```math
  E_{\mathrm{dmi}} = \int_V D \vec{m} \cdot (\nabla \times \vec{m}) \, \mathrm{d}V
  ```

- Anisotropy

```math
E_\mathrm{anis} = \int_{V} K_{u} [ 1 - (\vec{m} \cdot \hat{u})^2 ]\, dV
```

## LLG equation

The LLG equation is written as

```math
\frac{\partial \vec{m}}{\partial t} = - \gamma \vec{m} \times \vec{H} + \alpha \vec{m} \times  \frac{\partial \vec{m}}{\partial t}
```

and the corresponding LL form is given by

```math
(1+\alpha^2)\frac{\partial \vec{m}}{\partial t} = - \gamma \vec{m} \times \vec{H} - \alpha \gamma \vec{m} \times (\vec{m} \times \vec{H})
```

## LLG equation with extensions

For the driver `LLG_STT_CPP` the implemented equations is

```math
\frac{\partial \vec{m}}{\partial t} = - \gamma \vec{m} \times \vec{H} + \alpha \vec{m} \times  \frac{\partial \vec{m}}{\partial t}
+ (\vec{u} \cdot \nabla) \vec{m} - \beta [\vec{m}\times (\vec{u} \cdot \nabla)\vec{m}] - a_J \vec{m} \times (\vec{m} \times \vec{p})
 - \eta a_J \vec{m} \times \vec{p}
```

The simulation related to spin transfer torques (in-plane and current-perpendicular-to-plane) and the spin orbit torques can use the `LLG_STT_CPP` driver.

## Monte Carlo Simulation

For triangular mesh (2D), the system energy reads

```math
H= \sum_{\langle i, j\rangle}  \vec{D}_{i j} \cdot\left(\vec{S}_{i} \times \vec{S}_{j}\right)
-J \sum_{\langle i, j\rangle} \vec{S}_{i} \cdot \vec{S}_{j}- \lambda \sum_{\langle i, j\rangle} S_{i}^{z} S_{j}^{z}
-K \sum_{i}\left(S_{i}^{z}\right)^{2}
```

where

```math
\vec{D}_{i j} = D \hat{r}_{ij} \times \hat{z} + D_z^j \hat{z}
```
