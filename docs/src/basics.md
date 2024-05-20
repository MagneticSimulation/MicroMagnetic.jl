# Basics

## Meshes


MicroMag uses finite difference methods to discretize the micromagnetic energies. In MicroMag, the discretized grid information 
is stored in [`FDMesh`](@ref). Therefore, before starting the simulation, we need to create a mesh.

```julia
mesh = FDMesh(;dx=1e-9, dy=1e-9, dz=1e-9, nx=1, ny=1, nz=1)
```

In fact, [`FDMesh`](@ref) is used in our micromagnetic simulations, while for atomic models, we can use [`CubicMesh`](@ref), [`TriangularMesh`](@ref), as well as [`CylindricalTubeMesh`](@ref), etc.

```@raw html
<div class="mermaid">
graph LR;
    Mesh --> FDMesh
    Mesh --> AtomisticMesh
    AtomisticMesh --> TriangularMesh
    AtomisticMesh --> SquareMesh
    AtomisticMesh --> CubicMesh
    AtomisticMesh --> CylindricalTubeMesh
</div>
```

## Sim
In MicroMag, we can create different Sim objects based on the computational system or problem type. MicroMag defines four types of Sims:

```@raw html
<div class="mermaid">
graph LR
   AbstractSim --> MicroSim
   AbstractSim --> AtomisticSim
   AbstractSim --> NEB
   AbstractSim --> MonteCarlo
</div>
```

For [MicroSim](@ref MicroMag.MicroSim) and [AtomisticSim](@ref MicroMag.AtomisticSim), we recommend using the [`create_sim`](@ref) function to create them, as the [`create_sim`](@ref) function 
can specify some parameters while creating Sim. Of course, these parameters can also be specified later. 

```julia
sim = create_sim(mesh)
```
Note: all simulation data can be obtained through sim, especially, we can obtain the magnetization distribution state of the system through `sim.spin` at any time.

!!! note
    By default, the magnetization is stored in a 1D array with the form ``[m_{1,x}, m_{1, y}, m_{1, z}, ..., m_{n,x}, m_{n, y}, m_{n, z}]``, which can be reshaped into a 4D array
    ```julia
    m = reshape(sim.spin, 3, nx, ny, nz)
    mx = m[1, :, :, :]
    my = m[2, :, :, :]
    mz = m[3, :, :, :]
    ```


## Functions

In MicroMag, all parameters can be set using functions. For example, we can use the [`set_Ms`](@ref) function to set the saturation magnetization of the system. Of course, Ms should be a scalar for the same material, and we can set it like this:
```julia
set_Ms(sim, 8.6e5)
```
Additionally, we can set it with a function, like this:
```
function circular_Ms(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8.6e5
    end
    return 0.0
end
set_Ms(sim, circular_Ms)
```
Note that the Mesh we create is actually a regular cuboid, but in reality, the shape of the sample is not necessarily a cuboid. At this time, we define a round disk, 
where its Ms is 0 outside the disk. In this way, we can define the shape of the simulation system. Please note that in MicroMag, almost all setting functions can 
accept a function as input. This cell-based approach maximizes flexibility, allowing for defining shapes, defining multiple materials, etc.

## Shapes

### Basic Shapes

In addition to using functions to define shapes, for some regular shapes and their combinations, we can use basic shapes and boolean operations defined in MicroMag to achieve this. MicroMag supports Plane, Cylinder, Sphere, Box, and Torus, etc., as basic shapes.

!!! note 
    | **operator** | **Boolean operation** |
    | :----------: | :-------------------: |
    | +            | Union                 |
    | -            | Difference            |
    | *            | Intersection          |

Example:
```julia
using MicroMag

mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=50)

p1 = Plane(point=(40e-9,0,0), normal=(1, 0, 0))
save_vtk(mesh, p1, "shape1")

c1 = Cylinder(radius=30e-9, normal=(0.3,0,1))
save_vtk(mesh, c1, "shape2")

s1 = Sphere(radius = 30e-9, center=(50e-9, 0, 0))
save_vtk(mesh, s1, "shape3")

b1 = Box(sides = (110e-9, 50e-9, Inf), theta=pi/4)
save_vtk(mesh, b1, "shape4")

t1 = Torus(R = 60e-9, r=20e-9)
save_vtk(mesh, t1, "shape5")

t2 = t1 - b1 
save_vtk(mesh, t2, "shape6")

t3 = t1 * b1 
save_vtk(mesh, t3, "shape7")

t4 = t1 - p1 + (s1 * p1)
save_vtk(mesh, t4, "shape8")
```
The saved vts files can be visualized using programs such as Paraview, as shown below:

![shapes](./figures/shapes.png)

The created shapes can be used to set parameters, such as
```julia
set_Ms(sim::AbstractSim, geo::Shape, Ms::Number)
```
### Custom Shapes

We can also define custom shapes using the `create_shape` function. Custom shapes can also be combined with basic shapes using boolean operations.

## Energy Terms

After creating Sim, we can call functions to add energy terms that need to be considered in the simulation. For example, add_zeeman, add_exch, add_dmi, add_demag 
correspondingly add Zeeman energy, exchange interaction energy, DMI, and demagnetization energy.

Note
```julia
sim = create_sim(mesh, Ms=8e5, A=1.3e-11)
```
and
```julia
sim = create_sim(mesh)
set_Ms(sim, Ms=8e5)

ex = add_exch(sim, A=1.3e-11)
```
are equivalent. The advantage of the latter is that when we need exchange interaction data, we can directly access it through  `ex`.

MicroMag implements energy terms 

```@raw html
<div class="mermaid">
graph TD;
    MicroEnergy --> Exchange
    Exchange --> UniformExchange
    Exchange --> SpatialExchange
    Exchange --> ExchangeRKKY
    MicroEnergy --> BulkDMI
    MicroEnergy --> SpatialBulkDMI
    MicroEnergy --> Zeeman
    MicroEnergy --> Anisotropy
    MicroEnergy --> CubicAnisotropy
    MicroEnergy --> StochasticField    
</div>
```


## Driver

```@raw html
<div class="mermaid">
graph LR;
    Driver --> LLG
    Driver --> LLG_STT
    Driver --> EnergyMinimization
</div>
```



## Periodic Boundary conditions