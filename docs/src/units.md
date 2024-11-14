# Units

MicroMagnetic.jl uses **SI units** by default in simulations.

## Micromagnetic Model

In the Micromagnetic model, all parameters are expressed using SI units. Below is a summary of the physical quantities and their corresponding units.

| Physical Quantity           | Units   | Example Usage                                                 |
|:---------------------------:|:-------:|:-------------------------------------------------------------:|
| **Length**                 | m       | `mesh = FDMesh(dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1)` |
| **Time**                   | s       | `run_sim(sim; steps=100, dt=1e-11)`                           |
| **Saturation Magnetization** | A/m     | `set_Ms(sim, 8.6e5)`                                          |
| **Exchange Constant**      | J/m     | `add_exch(sim, 1.3e-11)`                                      |
| **DMI Constant**           | J/m²    | `add_dmi(sim, 1e-3)`                                          |
| **Anisotropy Constant**    | J/m³    | `add_anis(sim, 1e5; axis=(0,0,1))`                            |
| **External Magnetic Field** | A/m     | `add_zeeman(sim, (-24.6mT, 4.3mT, 0))`                        |
| **Gyromagnetic Ratio**     | m/(A·s) | `sim.driver.gamma = 2.21e5`                                   |
| **Temperature**            | K       | `add_thermal_noise(sim, 100.0)`                               |

## Atomistic Model

The Atomistic model also use SI units where applicable but often expressed in terms of energy or magnetic moments. 

| Physical Quantity         | Units    | Example Usage                                                 |
|:-------------------------:|:--------:|:-------------------------------------------------------------:|
| **Length**               | m        | `mesh = CubicMesh(dx=1e-9, dy=1e-9, dz=1e-9, nx=10, ny=10, nz=1)` |
| **Time**                 | s        | `run_sim(sim; steps=100, dt=1e-11)`                           |
| **Magnetic Moment**      | A·m²     | `set_mu_s(sim, mu_s_1)`                                       |
| **Exchange Interaction** | J        | `add_exch(sim, 50 * k_B)`                                     |
| **DMI Constant**         | J        | `add_dmi(sim, 5 * k_B)`                                       |
| **Anisotropy Energy**    | J        | `add_anis(sim, 5 * k_B; axis=(0,0,1))`                        |
| **External Magnetic Field** | T     | `add_zeeman(sim, (0, 0, 0.1))`                                |
| **Gyromagnetic Ratio**   | rad/(T·s) | `sim.driver.gamma = 1.76e11`                                  |
| **Temperature**          | K        | `add_thermal_noise(sim, 100.0)`                               |

## Dimensionless Units

In certain cases, simulations can be simplified using **dimensionless units**, particularly useful in theoretical studies or scaling analyses. 

### Example of Dimensionless Setup

| Parameter               | Value          |
|:-----------------------:|:--------------:|
| **Lattice Constant**    | $a = 1$        |
| **Spin Length**         | $S = 1$       |
| **Gyromagnetic Ratio**  | $ \gamma = 1 $ |
| **Magnetic Moment**     | $ \mu_s = 1 $ |
| **Exchange Constant**   | $ J = 1 $    |
| **DMI Strength**        | $ D/J = 0.09 $ |
| **Magnetic Field**      | $ H_0 = 0.00729 $ |
| **Gilbert Damping**     | $ \alpha = 0.04 $ |


Below is an example of using dimensionless units for setting up a simulation in MicroMagnetic.jl:

```julia
using MicroMagnetic

# Define a cubic mesh with normalized lattice units
mesh = CubicMesh(; nx=200, ny=200, nz=1, dx=1, dy=1, dz=1)

# Initialize the simulation with the Landau-Lifshitz-Gilbert (LLG) driver
sim = Sim(mesh; driver="LLG", name="test")
sim.driver.gamma = 1
sim.driver.alpha = 0.04

# Set magnetic moment
set_mu_s(sim, 1)

# Initialize magnetization direction
init_m0(sim, (1, 0.2, 0))

# Add exchange and DMI interactions
add_exch(sim, 1; name="exch")
add_dmi(sim, 0.09; name="dmi")

# Apply a dimensionless external magnetic field
add_zeeman(sim, (0, 0, 0.00729))
```