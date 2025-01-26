# Contributing

### Step-by-Step Guide

#### 1. **Fork the Repository**
Forking creates a personal copy of the MicroMagnetic.jl repository under your GitHub account.

- Go to the **MicroMagnetic.jl** repository: [MicroMagnetic.jl](https://github.com/magneticsimulation/MicroMagnetic.jl).
- In the top-right corner, click the **Fork** button to create a copy of the repository in your own GitHub account.

#### 2. **Clone Your Fork**
After forking the repository, you need to clone your fork to your local machine to start working on it.

- Open your terminal (or Git Bash on Windows).
- Clone the repository using the following command (replace `yourusername` with your GitHub username):
  ```bash
  git clone https://github.com/yourusername/MicroMagnetic.jl
  cd MicroMagnetic.jl
  ```

#### 3. **(Optional) Set the Original Repository as Upstream**
To keep your fork in sync with the original repository, you can add the original `MicroMagnetic.jl` repository as a remote named `upstream`. This step is optional, but it helps you easily pull in updates from the main repository.

- Run the following commands to set it up:
  ```bash
  git remote add upstream https://github.com/magneticsimulation/MicroMagnetic.jl
  git fetch upstream
  ```

If you don't set up `upstream`, you can still contribute, but you’ll need to manually manage any updates from the original repository.

#### 4. **Create a New Branch**
Before making any changes, create a new branch to isolate your modifications from the main branch. It's good practice to create descriptive branch names, like `fix-simulation-bug` or `add-feature-x`.

- Create and switch to a new branch:
  ```bash
  git checkout -b my-new-feature
  ```

#### 5. **Activate Local Development in Julia**
To ensure that Julia uses your local version of the package rather than the registered one, activate development mode in Julia. Assuming you have cloned the repository to `/path/to/MicroMagnetic.jl`, use the following command:

```julia
  using Pkg
  Pkg.develop(path="/path/to/MicroMagnetic.jl")
```

This allows you to work on the local version of the package.

#### 6. **Make Your Changes**
Now that you are in your new branch, you can start making changes to the code. For example, you can add new functionality, fix bugs, or update documentation.

- Use your favorite text editor or IDE to modify the source code in the **MicroMagnetic.jl** repository.

#### 7. **Test Your Changes**
Before pushing your changes, ensure they work correctly by running tests. The tests are usually located in the `test` folder.

- To run tests:
  ```bash
  julia --project test/runtests.jl
  ```

Make sure all tests pass before proceeding. If any tests fail, debug and fix your code.

#### 8. **Commit Your Changes**
Once you're happy with the changes, you need to commit them. Make sure to write a clear and concise commit message describing what you’ve done.

- Stage your changes:
  ```bash
  git add .
  ```
- Commit with a message:
  ```bash
  git commit -m "Add feature X to improve simulation accuracy"
  ```

#### 9. **Push Your Changes to Your Fork**
Now that you’ve committed your changes locally, push them to your forked repository on GitHub.

- Push your branch:
  ```bash
  git push origin my-new-feature
  ```

#### 10. **Open a Pull Request**
Once your changes are pushed to your fork on GitHub, you’re ready to open a pull request (PR) to merge your changes back into the main repository.

- Go to your forked repository on GitHub.
- You’ll see a button prompting you to open a pull request. Click on it.
- In the pull request, provide a detailed description of the changes you made, why they are needed, and any relevant information.
- Submit the pull request.

#### 11. **(Optional) Keep Your Fork Up-to-Date**
As the original `MicroMagnetic.jl` repository evolves, you’ll want to keep your fork in sync with the latest changes. If you followed step 3 to set `upstream`, you can do this easily.

- Fetch the latest updates from the original repository:
  ```bash
  git fetch upstream
  ```
- Merge the changes into your local branch:
  ```bash
  git merge upstream/main
  ```

### An Example: Adding Hexagonal Anisotropy

All changes can be found at [this commit](https://github.com/MagneticSimulation/MicroMagnetic.jl/commit/6a9d7ec7f01cacfaa7a908dc286f58a5bcae910b).

#### 1. **Find the Energy Expression**
The energy density of hexagonal anisotropy is expressed as:

```math
E = K_1 \sin^2 \theta + K_2 \sin^4 \theta + K_3 \sin^6 \theta \cos 6\phi
```

Here, $\theta$ is the angle between the magnetization vector $\mathbf{m}$ and the $c$-axis ($z$-axis).
$\phi$ is the angle of the projection of $\mathbf{m}$ on the hexagonal plane, measured with respect to the $x$-axis.

#### 2. **Compute the Effective Field**
Using the identity $\cos 6x = -\sin^6 x + 15 \cos^2 x \sin^4 x - 15 \cos^4 x \sin^2 x + \cos^6 x$, the energy density can be rewritten in terms of $m_x$, $m_y$, and $m_z$ as:

```math
E = K_1 (1 - m_z^2) + K_2 (1 - m_z^2)^2 + K_3 \left( m_x^6 - 15m_x^4 m_y^2 + 15m_x^2 m_y^4 - m_y^6 \right)
```

The corresponding effective field is:

```math
\mathbf{H}_\mathrm{eff} = -\frac{6K_3}{\mu_0 M_s} \left( m_x^5 - 10m_x^3m_y^2 + 5m_xm_y^4 \right) \mathbf{e}_x 
 -\frac{6K_3}{\mu_0 M_s} \left( -6m_x^4m_y + 10m_x^2m_y^3 - m_y^5 \right) \mathbf{e}_y 
+ \frac{2m_z}{\mu_0 M_s} \left[ K_1 + 2K_2(1 - m_z^2) \right] \mathbf{e}_z
```

#### 3. **Define the Struct in `src/head.jl`**
```julia
mutable struct HexagonalAnisotropy{T<:AbstractFloat} <: MicroEnergy
    K1::T
    K2::T 
    K3::T
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end
```

#### 4. **Add the Interface in `src/micro/add_field.jl`**
Add the following interface and export it. Including documentation is also recommended.

```julia
function add_hex_anis(sim::AbstractSim, K1=0, K2=0, K3=0, name="hex")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    T = Float[]
    anis = HexagonalAnisotropy(T(K1), T(K2), T(K3), field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return anis
end
```

#### 5. **Implement the Kernel in `src/micro/kernels.jl`**
```julia
@kernel function hexagonal_anisotropy_kernel!(@Const(m), h, energy, K1::T, K2::T, K3::T, 
                                            @Const(mu0_Ms), volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)

    @inbounds Ms_local = mu0_Ms[id]

    if Ms_local == 0.0
        @inbounds energy[id] = 0
        @inbounds h[j + 1] = 0
        @inbounds h[j + 2] = 0
        @inbounds h[j + 3] = 0
    else
        Ms_inv::T = 1.0 / Ms_local
        @inbounds mx = m[j + 1]
        @inbounds my = m[j + 2]
        @inbounds mz = m[j + 3]
        @inbounds h[j + 1] = -6*K3*Ms_inv*(mx^5-10*mx^3*my^2+5*mx*my^4)
        @inbounds h[j + 2] = -6*K3*Ms_inv*(-5*mx^4*my+10*mx^2*my^3-my^5)
        @inbounds h[j + 3] = 2*mz*Ms_inv*(K1 + 2*K2*(1-mz*mz))
        @inbounds energy[id] = (K1*(1-mz*mz) + K2*(1-mz*mz)^2 + K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)) * volume
    end
end
```

#### 6. **Implement `effective_field` Function in `src/micro/field.jl`**
```julia
function effective_field(anis::HexagonalAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    hexagonal_anisotropy_kernel!(default_backend[])(spin, anis.field, anis.energy, anis.K1, 
                                            anis.K2, anis.K3, sim.mu0_Ms, volume; ndrange=N)

    return nothing
end
```

#### 7. **Add Unit Tests (`test/test_anis.jl`)**
In the `test_hex_anis` function, compute the effective field for a given magnetization and compare it with results obtained using automatic differentiation.

```julia
function hexagonal_energy(m, K1, K2, K3)
    mx, my, mz = m[1], m[2], m[3]
    return K1*(1-mz*mz) + K2*(1-mz*mz)^2 + K3*(mx^6-15*mx^4*my^2+15*mx^2*my^4-my^6)
end

function test_hex_anis()
    mesh = FDMesh(; nx=10, ny=1, nz=1)

    sim = Sim(mesh)
    Ms = 8.6e5
    set_Ms(sim, Ms)
    m0 = (0.7, -0.4, 1.2)
    init_m0(sim, m0; norm=false)

    K1, K2, K3 = 1.23e2, 3.7e3, 6.9e2
    anis = add_hex_anis(sim, K1=K1, K2=K2, K3=K3)

    MicroMagnetic.effective_field(sim, sim.spin, 0.0)

    field = Array(anis.field)
    energy = Array(anis.energy)

    gd = Enzyme.gradient(Forward, hexagonal_energy, m0, Const(K1), Const(K2), Const(K3))
    expected = - collect(gd[1]) ./ (MicroMagnetic.mu_0*Ms)

    @test isapprox(field[1:3], expected)
    @test isapprox(energy[1]*1e27, hexagonal_energy(m0, K1, K2, K3), rtol=1e-5)
end
```

#### 8. **Add Documentation**
It's recommended to document the implementation for better clarity and user understanding.


### Summary of Key Commands
Here’s a quick reference for the key Git commands used in this guide:

- Fork: (via GitHub UI)
- Clone: `git clone https://github.com/yourusername/MicroMagnetic.jl`
- (Optional) Set Upstream: `git remote add upstream https://github.com/magneticsimulation/MicroMagnetic.jl`
- New Branch: `git checkout -b my-new-feature`
- Add Changes: `git add .`
- Commit: `git commit -m "message"`
- Push: `git push origin my-new-feature`
- (Optional) Fetch Updates: `git fetch upstream`
- (Optional) Merge Updates: `git merge upstream/main`
- Activate Local Development: `Pkg.develop(path="/path/to/MicroMagnetic.jl")`

Feel free to explore, experiment, and contribute to making **MicroMagnetic.jl** even better!