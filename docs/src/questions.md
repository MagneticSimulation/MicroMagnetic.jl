# Questions

### How the magnetization is stored in the array `spin`?

In JuMag, the magnetization is stored in a 1D array with the form ``[m_{1,x}, m_{1, y}, m_{1, z}, ..., m_{n,x}, m_{n, y}, m_{n, z}]``

### How to get the global index of magnetization at site ``(i,j,k)`` ?

In JuMag, the global index can be obtained using the following function

```julia
function index(i::Int64, j::Int64, k::Int64, nx::Int64, ny::Int64, nz::Int64)
    if i < 1 || j < 1 || k < 1 || k > nz || j > ny || i > nx
        return -1
    end
    return (k-1) * nx*ny + (j-1) * nx + i
end
```

### How to get the effective field at site ``(i,j,k)``?

The effective field is stored in the same form of the magnetization, so the effective at site `(i,j,k)`
can be extracted using

```julia
  id = index(i,j,k, nx, ny, nz)
  fx = sim.field[3*id-2]
  fy = sim.field[3*id-1]
  fz = sim.field[3*id]
```

Alternatively, we can use reshape function

```julia
  f = reshape(sim.field, 3, nx, ny, nz)
  fx,fy,fz = f[:, i,j,k]
```

### How to run JuMag in a cluster node with multiple GPUs?

In the long run, JuMag will support multiple GPUs. However, here we discuss the common scenario that one needs to change
some parameters systematically, such as the geometry size, or cell size, or external field, or charge current as a driving force.
In these cases, the simulation can be run in parallel since they are independent tasks. Here is a demo that using 4 GPUs to relax
the system with different sizes. 

```julia
@everywhere using JuMag
@everywhere using Printf
@everywhere using CUDAnative

@everywhere function relax_system(Nx, gpu_id)

  device!(gpu_id)  #using the gpu with gpu_id

  mesh =  FDMeshGPU(nx=Nx, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9)

  name = @sprintf("Nx_%d", Nx)
  sim = Sim(mesh, name=name, driver="SD")
  set_Ms(sim, 8.0e5)
  sim.driver.min_tau = 1e-10

  add_exch(sim, 1.3e-11)
  add_demag(sim)

  init_m0(sim, (1, 0.25, 0.1))

  relax(sim, maxsteps=5000, stopping_torque=1.0)

  return nothing
end

Nx = [i for i in 200:20:400]
for i = 1:4:length(Nx)
    r1 = remotecall(relax_system, 1, Nx[i], 0)
    r2 = remotecall(relax_system, 2, Nx[i+1], 1)
    r3 = remotecall(relax_system, 3, Nx[i+2], 2)
    r4 = remotecall(relax_system, 4, Nx[i+3], 3)
    fetch(r1)
    fetch(r2)
    fetch(r3)
    fetch(r4)
end
```

Here is the corresponding Slurm script:
```
#!/bin/bash
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks=1            # Number of tasks?
#SBATCH --cpus-per-task=4     # Number of CPUs per task
#SBATCH --gres=gpu:4          # Number of GPUs
srun julia -p 4 main.jl
```
