export add_zeeman, update_zeeman, add_anis, update_anis, add_cubic_anis, add_hex_anis,
       add_exch, add_dmi, add_sahe_torque, add_demag, add_dmi_int, add_exch_int,
       add_thermal_noise, add_sot, add_stt, add_torque

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Add a static Zeeman energy to the simulation.
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    n_total = sim.n_total

    T = Float[]
    field = zeros(T, 3 * n_total)

    init_vector!(field, sim.mesh, H0)

    field_kb = kernel_array(field)
    energy_kb = create_zeros(n_total)
    if isa(H0, Tuple)
        zeeman = Zeeman(H0, field_kb, energy_kb, name)
    else
        zeeman = Zeeman((0, 0, 0), field_kb, energy_kb, name)
    end
    push!(sim.interactions, zeeman)

    if sim.save_data
        id = length(sim.interactions)
        if isa(H0, Tuple) && length(H0) == 3
            field_item = SaverItem((string(name, "_Hx"), string(name, "_Hy"),
                                    string(name, "_Hz")), ("<A/m>", "<A/m>", "<A/m>"),
                                   o::AbstractSim -> o.interactions[id].H0)
            push!(sim.saver.items, field_item)
        end
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end

    @info "Static Zeeman has been added."

    return zeeman
end

"""
    update_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Set the Zeeman field to H0 where H0 is TupleOrArrayOrFunction according to its name. For example,

```julia
   add_zeeman(sim, (0,0,0), name="my_H")  #create a zeeman energy with field (0,0,0) A/m
   update_zeeman(sim, (0,0,1e5), name="my_H")  #change the field to (0,0,1e5) A/m
```
"""
function update_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    n_total = sim.n_total
    T = Float[]
    field = zeros(T, 3 * n_total)
    init_vector!(field, sim.mesh, H0)

    for i in sim.interactions
        if i.name == name
            copyto!(i.field, field)
            if isa(H0, Tuple)
                i.H0 = H0
            end
            return nothing
        end
    end
    return nothing
end

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")

Add a time-varying Zeeman interaction to the simulation system.

# Arguments
- `sim::AbstractSim`: The simulation object to which the Zeeman interaction will be added.
- `H0::TupleOrArrayOrFunction`: The spatial field configuration, which can be:
  - A tuple of 3 numbers representing a uniform field in x, y, z directions
  - An array representing a spatially varying field
  - A function that returns the field at specific spatial coordinates
- `ft::Function`: The time-dependent function that modulates the field amplitude.
- `name::String`: (Optional) A name identifier for this interaction. Defaults to "timezeeman".

# Time Function Specification
The time function `ft` should accept a time argument `t` and return either:
- **Scalar return**: A single value that will be applied uniformly to all three field components
- **Tuple return**: A 3-element tuple `(fx, fy, fz)` specifying individual scaling factors for x, y, z components

The function is evaluated at each time step to modulate the applied field.

# Examples

## Example 1: Uniform rotating field with vector time function
```julia
function time_fun(t)
    w = 2*pi*2.0e9  # 2 GHz frequency
    return (sin(w*t), cos(w*t), 0)  # Rotating in xy-plane
end

function spatial_H(i, j, k, dx, dy, dz)
    H = 1e3  # 1000 A/m base field strength
    if i <= 2
        return (H, H, 0)  # Apply field only in first two cells
    end
    return (0, 0, 0)      # Zero field elsewhere
end

add_zeeman(sim, spatial_H, time_fun)
```

## Example 2: Uniformly oscillating field with scalar time function
```julia
function time_fun_scalar(t)
    w = 2*pi*1.0e9  # 1 GHz frequency
    return sin(w*t)  # Same modulation for all components
end

# Apply uniform field in z-direction
add_zeeman(sim, (0, 0, 1e3), time_fun_scalar)
```

## Example 3: Complex modulated field
```julia
function complex_time_fun(t)
    w1 = 2*pi*1.0e9
    w2 = 2*pi*0.5e9
    # Different modulation for each component
    fx = sin(w1*t) * cos(w2*t)
    fy = 0.5 * (1 + sin(w1*t))
    fz = exp(-t/1e-9)  # Exponential decay in z-component
    return (fx, fy, fz)
end

# Spatially uniform initial field
add_zeeman(sim, (1e3, 2e3, 0.5e3), complex_time_fun, name="modulated_zeeman")
```

# Field Initialization
The `H0` parameter defines the spatial distribution of the field:
- For uniform fields: Use a 3-element tuple `(Hx, Hy, Hz)`
- For spatially varying fields: Provide a function with signature `(i, j, k, dx, dy, dz) -> (Hx, Hy, Hz)`
  where `i, j, k` are grid indices and `dx, dy, dz` are spatial steps

# Energy Calculation
The Zeeman energy density is calculated as:
E = -μ₀ * M ⋅ H
where M is the magnetization and H is the applied field.
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function;
                    name="timezeeman")
    n_total = sim.n_total
    T = Float[]

    init_field = zeros(T, 3 * n_total)
    init_vector!(init_field, sim.mesh, H0)

    field_kb = kernel_array(init_field)
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    f0 = ft(0)
    ft_length = f0 isa Tuple ? length(f0) : 1
    if ft_length != 3 && ft_length != 1
        error("The time_fun function should return either 1 value (scalar) or 3 values (as a tuple).")
    end

    is_scalar = ft_length == 1 ? true : false
    zee = TimeZeeman(T(0), T(0), T(0), ft, is_scalar, field_kb, field, energy, name)
    push!(sim.interactions, zee)

    if sim.save_data
        unit = isa(sim, MicroSim) ? "<A/m>" : "<Tesla>"
        id = length(sim.interactions)
        if isa(H0, Tuple) && length(H0) == 3
            field_item = SaverItem((string(name, "_Hx"), string(name, "_Hy"),
                                    string(name, "_Hz")), (unit, unit, unit),
                                   o::AbstractSim -> (H0[1]*zee.time_fx, H0[2]*zee.time_fy,
                                                      H0[3]*zee.time_fz))
            push!(sim.saver.items, field_item)
        end
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return zee
end

@doc raw"""
    add_exch(sim::MicroSim, A::NumberOrTupleOrArrayOrFunction; name="exch")

Add exchange energy to the system. The exchange energy is definded as 
```math
  E_\mathrm{ex} = \int_{V} A (\nabla \mathbf{m})^2 \mathrm{d}V
```

# Examples:

```julia
    add_exch(sim, 1e-11)
```

or 

```julia
    add_exch(sim, (2e-12,5e-12,0))
```

or

```julia
    function spatial_A(i,j,k,dx,dy,dz)
        if i<10
            return 1e-11
        else
            return 2e-11
        end
    end
    add_exch(sim, spatial_A)
```
"""
function add_exch(sim::MicroSim, A::NumberOrTupleOrArrayOrFunction; name="exch")
    n_total = sim.n_total
    T = Float[]
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    exch = nothing
    if isa(A, Number)
        exch = UniformExchange(T(A), T(A), T(A), field, energy, name)
    elseif isa(A, Tuple) && length(A) == 3
        exch = UniformExchange(T(A[1]), T(A[2]), T(A[3]), field, energy, name)
    else
        Spatial_A = zeros(T, sim.n_total)
        init_scalar!(Spatial_A, sim.mesh, A)

        A_kb = kernel_array(Spatial_A)

        exch = SpatialExchange(A_kb, field, energy, name)
    end

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    @info "Exchange has been added."
    return exch
end

@doc raw"""
    add_dmi(sim::MicroSim, D::NumberOrTupleOrArrayOrFunction; name="dmi", type="bulk")

Add DMI to the system. `type` could be "bulk", "interfacial" or "D2d". 

Examples:

```julia
   add_dmi(sim, 1e-3, type="interfacial")
```
or
```julia
   add_dmi(sim, 1e-3, type="D2d")
```

```julia
   add_dmi(sim, (1e-3, 1e-3, 0), type="bulk")
```

"""
function add_dmi(sim::MicroSim, D::NumberOrTupleOrArrayOrFunction; name="dmi", type="bulk")
    n_total = sim.n_total
    T = Float[]
    field = KernelAbstractions.zeros(default_backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(default_backend[], T, n_total)

    if type == "bulk"
        if isa(D, Number)
            dmi = BulkDMI(T(D), T(D), T(D), field, energy, name)
        elseif isa(D, Tuple) && length(D) == 3
            dmi = BulkDMI(T(D[1]), T(D[2]), T(D[3]), field, energy, name)
        else
            Spatial_D = zeros(T, sim.n_total)
            init_scalar!(Spatial_D, sim.mesh, D)

            D_kb = KernelAbstractions.zeros(default_backend[], T, n_total)
            copyto!(D_kb, Spatial_D)

            dmi = SpatialBulkDMI(D_kb, field, energy, name)
        end

        @info "Bulk DMI has been added."
    elseif type == "interfacial"
        Spatial_D = zeros(T, sim.n_total)
        init_scalar!(Spatial_D, sim.mesh, D)

        D_kb = KernelAbstractions.zeros(default_backend[], T, n_total)
        copyto!(D_kb, Spatial_D)

        dmi = InterfacialDMI(D_kb, field, energy, name)
        @info "Interfacial DMI has been added."

    elseif type == "D2d"
        if isa(D, Number)
            dmi = BulkDMI(T(-D), T(D), T(0), field, energy, name)
        else
            error("D2d only support uniform DMI!")
        end
    else
        error("Supported DMI type:", " interfacial", " bulk", " D2d")
    end

    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return dmi
end

""" 
    add_exch(sim::AbstractSim, geo::Shape, A::Number; name="exch")

Add exchange interaction within the Shape, or update corresponding A if other exch is added.
"""
function add_exch(sim::MicroSim, geo::Shape, A::Number; name="exch")
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_Shape(interaction.A, geo, A)
            return nothing
        end
    end
    n_total = sim.n_total
    field = zeros(Float64, 3 * n_total)
    energy = zeros(Float64, n_total)
    Spatial_A = zeros(Float64, n_total)
    update_scalar_Shape(Spatial_A, geo, A)
    if isa(sim, MicroSim)
        exch = Exchange(Spatial_A, field, energy, name)
    else
        exch = HeisenbergExchange(A, field, energy, name)
    end
    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.ite,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return exch
end

"""
    add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0, fft=true)

Add Demag to the system. `Nx`, `Ny` and `Nz` can be used to describe the macro boundary conditions which means that
the given mesh is repeated `2Nx+1`, `2Ny+1 and `2Nz+1` times in `x`, `y` and `z` direction, respectively.
"""
function add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0, fft=true)
    if fft && Float[] != AbstractFloat
        demag = init_demag(sim, Nx, Ny, Nz)
    else
        demag = init_direct_demag(sim, Nx, Ny, Nz)
    end
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return demag
end

"""
    add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis::TupleOrArrayOrFunction, name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
    E_\\mathrm{anis} =  K_{u} (1 - \\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction;
                  axis::TupleOrArrayOrFunction=(0, 0, 1), name="anis")
    n_total = sim.n_total
    T = Float[]
    Kus = zeros(T, n_total)
    init_scalar!(Kus, sim.mesh, Ku)

    Kus_kb = kernel_array(Kus)
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    if isa(axis, Tuple)
        lt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)
        naxis = (T(axis[1] / lt), T(axis[2] / lt), T(axis[3] / lt))
        anis = Anisotropy(Kus_kb, naxis, field, energy, name)
    else
        axes = zeros(T, 3 * n_total)
        init_vector!(axes, sim.mesh, axis)
        normalise(axes, n_total)

        axes_kb = kernel_array(axes)
        anis = SpatialAnisotropy(Kus_kb, axes_kb, field, energy, name)
    end

    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end

    @info "Uniaxial Anisotropy has been added."
    return anis
end

"""
    add_anis(sim::AbstractSim, geo::Shape, Ku::Number; axis=(0,0,1), name="anis")

Add Anisotropy within the Shape, or update corresponding Ku if other anis is added.
"""
#FIXME : fix this function
function add_anis(sim::AbstractSim, geo::Shape, Ku::Number; axis=(0, 0, 1), name="anis")
    lt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)

    naxis = (axis[1] / lt, axis[2] / lt, axis[3] / lt)
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.Ku, geo, Ku)
            return nothing
        end
    end
    n_total = sim.n_total
    Kus = zeros(Float64, n_total)
    update_scalar_geometry(Kus, geo, Ku)
    field = zeros(Float64, 3 * n_total)
    energy = zeros(Float64, n_total)
    anis = Anisotropy(Kus, naxis, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return anis
end

"""
    update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name = "anis")

update anisotropy constant Ku according to its name.

Example:
```julia
    mesh = FDMesh(nx=200, ny=200, nz=12, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh)
    add_anis(sim, 3e4, axis = (0,0,1), name="K1")
    add_anis(sim, 1e5, axis = (1,0,0), name="K2")
    update_anis(sim, 5e4, name="K2")  #update anisotropy K2
```
"""
function update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name="anis")  #FIXME : fix this function
    n_total = sim.n_total
    Kus = zeros(Float64, n_total)
    init_scalar!(Kus, sim.mesh, Ku)
    for i in sim.interactions
        if i.name == name
            i.Ku[:] = Kus[:]
            return nothing
        end
    end
    return nothing
end

@doc raw"""
    add_cubic_anis(sim::AbstractSim, Kc::NumberOrArrayOrFunction; axis1=(1,0,0), axis2=(0,1,0), name="cubic")

add a cubic anisotropy with default axis (1,0,0) , (0,1,0), and (0,0,1). The third axis is defined as axis3 = axis1 x axis2.

```math
  E_\mathrm{cubic} = -\int_{V} K_c (m_x^4 + m_y^4 + m_z^4) \, dV
```

# Example
```julia
    add_cubic_anis(sim, 1e3, (1, 1, 0), (1, -1, 0))
```
"""
function add_cubic_anis(sim::AbstractSim, Kc::NumberOrArrayOrFunction; axis1=(1, 0, 0),
                        axis2=(0, 1, 0), name="cubic")
    n_total = sim.n_total
    T = Float[]
    Kcs = zeros(T, n_total)
    init_scalar!(Kcs, sim.mesh, Kc)

    norm1 = sqrt(axis1[1]^2 + axis1[2]^2 + axis1[3]^2)
    norm2 = sqrt(axis2[1]^2 + axis2[2]^2 + axis2[3]^2)
    naxis1, naxis2 = axis1 ./ norm1, axis2 ./ norm2
    if abs.(sum(naxis1 .* naxis2)) > 1e-10
        @error("The axis1 and axis2 are not perpendicular to each other!")
        return nothing
    end
    naxis3 = cross_product(axis1, axis2)

    Kcs_kb = KernelAbstractions.zeros(default_backend[], T, n_total)
    field = KernelAbstractions.zeros(default_backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(default_backend[], T, n_total)

    copyto!(Kcs_kb, Kcs)

    anis = CubicAnisotropy(Kcs_kb, naxis1, naxis2, naxis3, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return anis
end

@doc raw"""
    add_hex_anis(sim::AbstractSim; K1=0, K2=0, K3=0, name="hex")

Add hexagonal anisotropy to a simulation. The energy density of the anisotropy is defined as:

```math
E = K_1 \sin^2 \theta + K_2 \sin^4 \theta + K_3 \sin^6 \theta \cos 6\phi
```

# Example
```julia
add_hex_anis(sim, K1=1e3, K2=0, K3=1e2)
```
"""
function add_hex_anis(sim::AbstractSim; K1=0, K2=0, K3=0, name="hex")
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

@doc raw"""
    add_thermal_noise(sim::AbstractSim, Temp::NumberOrArrayOrFunction; name="thermal", T0=0, scaling=t -> 1.0, k_B=k_B)

Add thermal noise fields to the simulation based on specified temperature profile. 

The effective temperature is defined differently based on the temperature specification:

### Case 1: Spatially Defined Temperature (Constant, Array, or 3-parameter Function)
For temperature specifications that depend only on space:
```math
T_{\text{eff}}(x,y,z,t) = \text{Temp}(x,y,z) \times \text{scaling}(t) + \text{T0}
```
- `Temp(x,y,z)` defines the spatial distribution
- `scaling(t)` provides time modulation
- `T0` adds a constant offset

### Case 2: Spatiotemporally Defined Temperature (4-parameter Function)
For functions with explicit space and time dependence `f(x,y,z,t)`:
```math
T_{\text{eff}}(x,y,z,t) = \text{Temp}(x,y,z,t)
```
In this case, `scaling` and `T0` parameters are ignored.

For micromagnetic model, the thermal noise is defined as

```math
\mathbf{b}^u = \eta \sqrt \frac{2 \alpha k_B T}{\mu_0 M_s \gamma \Delta V \Delta t}
```
and $\eta \sim \mathcal{N}(0,1)$ is a standard normal random variable, $\gamma=2.211\times 10^5$ m/(A·s) is the gyromagnetic ratio,
$\Delta V$ is the cell volume and $\Delta t$ is the time step.

For the atomistic model, the thermal noise is defined as

```math
\mathbf{b}^u = \eta \sqrt \frac{2 \alpha k_B T}{\gamma \mu_s \Delta t}.
```
where $\eta \sim \mathcal{N}(0,1)$ is a standard normal random variable, $\gamma=1.76\times 10^{11}$ rad/(T·s),
and $\mu_s$ is the effective magnetic moment.

## Arguments
- `sim::AbstractSim`: Simulation object to which thermal noise will be added
- `Temp::NumberOrArrayOrFunction`: Temperature specification. Can be:
  - Constant number (uniform temperature)
  - Array (spatially varying temperature)
  - 3-parameter function `f(x,y,z)` for spatial control
  - 4-parameter function `f(x,y,z,t)` for spatiotemporal control
- `name::String`: Identifier for the noise field (default: `"thermal"`)
- `T0::Number`: Base temperature offset (default: `0`)
- `scaling::Function`: Time-dependent scaling function `f(t)` (default: `t -> 1.0`)
- `k_B::Float64`: Boltzmann constant (default: global `k_B`)

## Examples

```julia
# 300 K everywhere
add_thermal_noise(sim, 300.0) 

# Constant temperature with time modulation
add_thermal_noise(sim, 100.0, scaling=t -> exp(-t/1e-9))

# Spatially varying temperature with time modulation
gaussian_profile = (x, y, z) -> 200 * exp(-(x^2 + y^2)/1e-18)
add_thermal_noise(sim, gaussian_profile, T0=50, scaling=t -> 0.5*sin(2π*t/1e-9) + 0.5)

# Fully spatiotemporal temperature (scaling and T0 ignored)
dynamic_temp = (x, y, z, t) -> 300 + 20*sin(2π*t/1e-9 + 0.1*x)*exp(-t/2e-9)
add_thermal_noise(sim, dynamic_temp)
```
"""
function add_thermal_noise(sim::AbstractSim, Temp::NumberOrArrayOrFunction; name="thermal",
                           T0=0, scaling=t -> 1.0, k_B=k_B)
    N = sim.n_total
    T = Float[]
    field = create_zeros(3 * N)
    energy = create_zeros(N)

    Spatial_T = create_zeros(N)
    eta = create_zeros(3 * N)

    if isa(Temp, Function) && methods(Temp)[1].nargs - 1 == 4 #four parameters (x,y,z,t)
        @info("Fully spatiotemporal temperature profile is provided, scaling and T0 are ignored.")
        spatiotemporal_mode = true
        scaling_fun = Temp
    else
        spatiotemporal_mode = false
        init_scalar!(Spatial_T, sim.mesh, Temp)
        scaling_fun = scaling
    end

    thermal = StochasticField(Spatial_T, T(T0), eta, field, energy, -1, name, k_B,
                              scaling_fun, average(Spatial_T), T(1), spatiotemporal_mode)
    push!(sim.interactions, thermal)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))

        push!(sim.saver.items,
              SaverItem(string("T_", name), "<K>",
                        o::AbstractSim -> o.interactions[id].mean_temperature *
                                          o.interactions[id].scaling_factor))
    end
    return thermal
end

@doc raw"""
    add_exch_int(sim::AbstractSim, J::NumberOrArrayOrFunction; k1=1, k2=-1, name="rkky")

Add an RKKY-type exchange for interlayers. The energy of RKKY-type exchange is defined as 

```math
E_\mathrm{rkky} =  - \int_\Gamma J_\mathrm{rkky} \mathbf{m}_{i} \cdot \mathbf{m}_{j} dA
```
where $\Gamma$ is the interface between two layers with magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$,
$J_\mathrm{rkky}$ is the coupling constant which is related to the spacer layer thickness. 

The effective field is given then as
```math
\mathbf{H}_i = \frac{1}{\mu_0 M_s}  \frac{J_\mathrm{rkky}}{\Delta_z} \mathbf{m}_{j} 
```
"""
function add_exch_int(sim::MicroSim, J::NumberOrArrayOrFunction; k1=1, k2=-1,
                      name="exch_int")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    mesh = sim.mesh

    mesh_tmp = FDMesh(nx=mesh.nx, ny=mesh.ny, nz=1, dx=mesh.dx, dy=mesh.dy, dz=mesh.dz)

    Js = create_zeros(mesh.nx*mesh.ny)
    init_scalar!(Js, mesh_tmp, J)

    if k2 == -1
        k2 = sim.mesh.nz
    end
    exch = InterlayerExchange(Js, Int32(k1), Int32(k2), field, energy, name)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return exch
end

@doc raw"""
    add_dmi_int(sim::MicroSimGPU, D::Tuple{Real, Real, Real}; k1=1, k2=-1, name="dmi_int")

Add an interlayer DMI to the system. The energy of interlayer DMI is defined as 

```math
E_\mathrm{dmi-int} =  \int_\Gamma \mathbf{D} \cdot \left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right) dA
```
where $\Gamma$ is the interface between two layers with magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$. 
$\mathbf{D}$ is the effective DMI vector. 

The effective field is given
```math
\mathbf{H}_i = \frac{1}{\mu_0 M_s \Delta_z}  \mathbf{D} \times \mathbf{m}_{j} 
```
"""
function add_dmi_int(sim::MicroSim, D::Tuple{Real,Real,Real}; k1=1, k2=-1, name="dmi_int")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    if k2 == -1
        k2 = sim.mesh.nz
    end

    T = Float[]
    dmi = InterlayerDMI(T(D[1]), T(D[2]), T(D[3]), Int32(k1), Int32(k2), field, energy,
                        name)

    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return dmi
end

@doc raw"""
    add_sahe_torque(sim::AbstractSim, sigma_s::TupleOrArrayOrFunction, sigma_sa::TupleOrArrayOrFunction; a1=(1,0,0), a2=(0,1,0), beta=0, name="she")

Add an effective field to represent the torque due to spin anomalous Hall effect, which is given by 

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t} 
- \mathbf{m} \times (\mathbf{m} \times \boldsymbol{\sigma}_{x y}^\mathrm{S H E}) -  \beta\mathbf{m} \times \boldsymbol{\sigma}_{x y}^\mathrm{S H E}
```

The equivalent effective field is
```math
\mathbf{H}_\mathrm{she} = (1/\gamma)(\beta \boldsymbol{\sigma}_{x y}^\mathrm{S H E} + \mathbf{m} \times \boldsymbol{\sigma}_{x y}^\mathrm{S H E})
```
where 
```math
\boldsymbol{\sigma}_{x y}^\mathrm{S H E}=\widetilde{\boldsymbol{\sigma}}_\mathrm{S} \times \hat{\mathbf{a}}_2
-\left({\mathbf{m}} \cdot\left(\hat{\mathbf{a}}_2 \times \hat{\mathbf{a}}_1\right)\right)^2\left(\widetilde{\boldsymbol{\sigma}}_\mathrm{SA} \times \hat{\mathbf{a}}_2\right)
```

### Example
```julia
    add_sahe_torque(sim, (0.1, 0.2, 0.3), (0.3, 0.4, 0.5), beta=0.1)
```

"""
function add_sahe_torque(sim::AbstractSim, sigma_s::TupleOrArrayOrFunction,
                         sigma_sa::TupleOrArrayOrFunction; a1=(1, 0, 0), a2=(0, 1, 0),
                         beta=0, name="she")
    n_total = sim.n_total
    mesh = sim.mesh

    T = Float[]
    field = create_zeros(T, 3 * n_total)
    sigma_s3 = zeros(T, 3 * n_total)
    sigma_sa3 = zeros(T, 3 * n_total)

    init_vector!(sigma_s3, sim.mesh, sigma_s)
    init_vector!(sigma_sa3, sim.mesh, sigma_sa)

    for i in 1:mesh.n_total
        j = 3*(i-1) + 1
        a = sigma_s3[j]
        b = sigma_s3[j+1]
        c = sigma_s3[j+2]
        sigma_s3[j] = cross_x(a, b, c, a2[1], a2[2], a2[3])
        sigma_s3[j+1] = cross_y(a, b, c, a2[1], a2[2], a2[3])
        sigma_s3[j+2] = cross_z(a, b, c, a2[1], a2[2], a2[3])

        a = sigma_sa3[j]
        b = sigma_sa3[j+1]
        c = sigma_sa3[j+2]
        sigma_sa3[j] = cross_x(a, b, c, a2[1], a2[2], a2[3])
        sigma_sa3[j+1] = cross_y(a, b, c, a2[1], a2[2], a2[3])
        sigma_sa3[j+2] = cross_z(a, b, c, a2[1], a2[2], a2[3])
    end

    c2 = cross_product(a2, a1)
    torque = SAHETorqueField(kernel_array(sigma_s3), c2, kernel_array(sigma_sa3), T(beta),
                             field, name)

    push!(sim.interactions, torque)

    return torque
end

@doc raw"""
    add_sot(sim::AbstractSim, aj::NumberOrArrayOrFunction, bj::Number, p::Tuple{Real,Real,Real}; name="sot")

Add damping-like and field-like torque to LLG equation, which is useful to model spin-orbit torques (SOT) or Slonczewski torque with constant $\epsilon(\theta)$.

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H}_\text{eff} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t}
- a_J \mathbf{m} \times (\mathbf{m} \times \mathbf{p}) -  b_J \mathbf{m} \times \mathbf{p}
```

The equivalent effective field is
```math
\mathbf{H}_\mathrm{stt} = (1/\gamma)(a_J \mathbf{m} \times \mathbf{p} +  b_J \mathbf{p})
```
"""
function add_sot(sim::AbstractSim, aj::NumberOrArrayOrFunction, bj::Number,
                 p::Tuple{Real,Real,Real}; name="sot")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)

    T = Float[]
    spatial_aj = zeros(T, sim.n_total)
    init_scalar!(spatial_aj, sim.mesh, aj)
    aj_zb = kernel_array(spatial_aj)

    torque = DFTorqueField(T(p[1]), T(p[2]), T(p[3]), aj_zb, T(bj), field, name)

    push!(sim.interactions, torque)

    return torque
end

@doc raw"""
    add_torque(sim::AbstractSim, fun::Function; name="torque")

Add a torque $\mathbf{T}$ to the Landau-Lifshitz-Gilbert (LLG) equation.

The modified LLG equation becomes:

```math
\frac{\partial \mathbf{m}}{\partial t} = - \gamma \mathbf{m} \times \mathbf{H}_\text{eff} + \alpha \mathbf{m} \times  \frac{\partial \mathbf{m}}{\partial t} + \mathbf{T}
```

where $\mathbf{T} = \mathbf{T}(\mathbf{m}, t)$ is the additional torque term.

# Arguments
- `sim::AbstractSim`: The simulation object to which the torque will be added.
- `fun::Function`: A function that defines the torque term. The function should have the signature `fun(y, m, t)` where:
  - `y`: The torque
  - `m`: The magnetization vector m
  - `t`: The current simulation time

# Examples
```julia
function torque_fun(y, m, t)
    y .= 0 
    y[1:3:end] .= m[1:3:end] .* sin(t) # x-component
    y[2:3:end] .= 0.2 * cos(t) # y-component
end
add_torque(sim, torque_fun, name="test")
```
"""
function add_torque(sim::AbstractSim, fun::Function; name="torque")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)

    torque = TorqueField(fun, field, name)

    push!(sim.interactions, torque)
    return torque
end

@doc raw"""
    add_stt(sim::AbstractSim; model::Symbol=:slonczewski, name="stt", kwargs...)

Add spin transfer torque (STT) to the Landau-Lifshitz-Gilbert equation using various theoretical models.

Spin transfer torque is a fundamental phenomenon in spintronics where spin-polarized electric currents exert torques on magnetic moments, 
enabling precise manipulation of magnetization dynamics in nanoscale devices.

# Supported Models

## Zhang-Li Model (`model=:zhang_li`)

The Zhang-Li model describes spin transfer torque for in-plane current flows in continuous magnetic materials. 
This model is particularly suitable for describing current-induced domain wall motion, magnetization dynamics 
in nanowires, and spatially non-uniform current distributions.

The modified LLG equation becomes:

```math
\frac{d \mathbf{m}}{dt} = -\gamma \mathbf{m} \times \mathbf{H}_\mathrm{eff} + \alpha \mathbf{m} \times \frac{d \mathbf{m}}{dt}
- b \mathbf{m} \times [\mathbf{m} \times (\mathbf{j} \cdot \nabla) \mathbf{m}] 
- \xi b \mathbf{m} \times (\mathbf{j} \cdot \nabla) \mathbf{m}
```

where $\mathbf{j}$ is the current density vector and $\xi$ is the non-adiabatic parameter. The coefficient $b$ is calculated as:
```math
b = \frac{P \mu_B}{e M_s (1 + \xi^2)}
```

### Parameters

- **Required:**
  - `J::TupleOrArrayOrFunction`: Current density specification (A/m²)
    - **Constant vector:** `Tuple` or `Vector` specifying uniform current flow
    - **Spatially varying:** `Function` that returns current density at position `(x,y,z)`

- **Coefficient specification (choose one approach):**
  - **Physical parameters approach:**
    - `P::Real`: Spin polarization (0 ≤ P ≤ 1)
    - `Ms::Real`: Saturation magnetization (A/m)
  - **Direct coefficient approach:**
    - `b::Real`: Pre-calculated Zhang-Li torque coefficient

- **Optional:**
  - `xi::Real=0.0`: Non-adiabatic parameter (default: 0.0)
  - `ft::Function`: Time-dependent function modulating current density amplitude

### Examples

```julia
# Basic usage with uniform current vector
add_stt(sim, model=:zhang_li, P=0.5, Ms=8e5, xi=0.05, J=(1e12, 0, 0))

# Using pre-calculated coefficient
add_stt(sim, model=:zhang_li, b=72.5, xi=0.1, J=(1, 0, 0))

# Time-dependent current pulse with spatial variation
function current_pulse(t)
    return t < 1e-9 ? 1.0 : 0.0  # 1 ns pulse
end

function spatial_current(i, j, k, dx, dy, dz)
    current =  k > 5 ? 1e12 : 2e11
    return (current, 0.0, 0.0)
end
add_stt(sim, model=:zhang_li, P=0.5, Ms=8e5, J=spatial_current, ft=current_pulse)
```

## Slonczewski Model (`model=:slonczewski`)

The Slonczewski model describes spin transfer torque for perpendicular current flows in magnetic tunnel junctions (MTJs) 
and spin valves. The extended LLG equation is:

```math
\frac{d\mathbf{m}}{dt} = -\gamma \mathbf{m} \times \mathbf{H}_{\text{eff}}
   + \alpha \left( \mathbf{m} \times \frac{d\mathbf{m}}{dt} \right)
   + \gamma \beta \epsilon \left[ \mathbf{m} \times (\mathbf{m}_p \times \mathbf{m}) \right]
   - \gamma \beta \xi \mathbf{m} \times \mathbf{m}_p
```
where $\mathbf{m}=\mathbf{M}/M_s$, $\gamma$ the gyromagnetic ratio, $\mathbf{m}_p$ is electron polarization direction, and $\xi$ 
is the secondary spin-torque parameter. The coefficients $\beta$ and $\epsilon$ are defined as follows:
```math
\beta =  \frac{\hbar}{\mu_0 e} \frac{J}{t M_s}, \qquad \epsilon = \frac{P \Lambda^2}{(\Lambda^2 + 1) + (\Lambda^2 - 1)(\mathbf{m} \cdot \mathbf{m}_p)}
```
where $e$ is electron charge (C), $J$ is current density (A/m²), $t$ is free layer thickness (m), $M_s$ is the saturation magnetization (A/m),
$P$ is the spin polarization, $\Lambda$ is the Slonczewski parameter.

### Parameters
- **Required:**
  - `J::NumberOrArrayOrFunction`: Current density magnitude (A/m²)
  - `tf::Real`: Free layer thickness (m)
  - `Ms::Real`: Saturation magnetization (A/m)
  - `p::Union{Tuple, Vector}`: Polarization direction vector
  - `P::Real`: Spin polarization (0 ≤ P ≤ 1)

- **Optional:**
  - `xi::Real=0.0`: Secondary spin-torque parameter (default: 0.0)
  - `Lambda::Real=2.0`: Slonczewski parameter (default: 2.0)
  - `ft::Function`: Time-dependent function modulating current density amplitude

### Examples

```julia
# Const current density
add_stt(sim, model=:slonczewski, J=1e12, tf=2e-9, Ms=8e5, p=(0,0,1), P=0.7)

# With secondary spin-torque parameter
add_stt(sim, model=:slonczewski, J=8e11, tf=1.5e-9, Ms=6e5, p=(1,0,0), P=0.5, xi=0.05)


# Time-dependent current pulse with Gaussian current profile
function current_pulse(t)
    return t < 1e-9 ? 1.0 : 0.0  # 1 ns pulse
end
function gaussian_current(i, j, k, dx, dy, dz)
    x, y = (i-100)*dx, (j-100)*dy
    sigma = 50e-9  # 50 nm width
    r_sq = x^2 + y^2
    return 1e12 * exp(-r_sq/(2*sigma^2))
end

add_stt(sim, model=:slonczewski, J=gaussian_current, tf=2e-9, Ms=8e5, p=(0,0,1), P=0.4, ft=current_pulse)
```
"""
function add_stt(sim::AbstractSim; model::Symbol=:slonczewski, kwargs...)
    params = Dict(kwargs)
    if !haskey(params, :J)
        throw(ArgumentError("Current density parameter `J` is required."))
    end

    if model == :zhang_li
        return add_zhang_li_torque(sim, "zhangli", params)
    elseif model == :slonczewski
        return add_slonczewski_torque(sim, "slonczewski", params)
    else
        supported_models = [:zhang_li, :slonczewski]
        throw(ArgumentError("Unsupported STT model: $model. Supported models are: $supported_models"))
    end
end

function add_zhang_li_torque(sim::AbstractSim, name, params::Dict)
    has_b = haskey(params, :b)
    has_P = haskey(params, :P) && haskey(params, :Ms)

    if has_b && has_P
        @warn("Parameters `(P, Ms)`` will be ignored since `b` is provided.")
    end

    xi = get(params, :xi, 0)
    b = if has_b
        params[:b]
    elseif has_P
        P = params[:P]
        Ms = params[:Ms]
        P*mu_B/(c_e*Ms)/(1+xi^2)
    else
        throw(ArgumentError("Zhang-Li model requires either `b` or `(P, Ms)` parameters"))
    end

    T = Float[]
    n_total = sim.n_total

    bJ_cpu = zeros(T, 3 * n_total)
    init_vector!(bJ_cpu, sim.mesh, params[:J])
    bJ_cpu .*= b
    bJ = kernel_array(bJ_cpu)

    field = create_zeros(3 * n_total)

    ft = haskey(params, :ft) ? params[:ft] : t -> 1.0

    torque = ZhangLiTorque(T(xi), bJ, field, ft, name)

    push!(sim.interactions, torque)

    return torque
end

function add_slonczewski_torque(sim::AbstractSim, name, params::Dict)
    for k in [:tf, :Ms, :p, :P]
        if !haskey(params, k)
            throw(ArgumentError("Parameter $k is not provided"))
        end
    end

    xi = get(params, :xi, 0)
    Lambda = get(params, :Lambda, 2.0)
    tf = params[:tf]
    Ms = params[:Ms]
    P = params[:P]
    p = params[:p]

    beta = h_bar/(mu_0*Ms*c_e*tf)

    T = Float[]
    n_total = sim.n_total

    J_cpu = zeros(T, n_total)
    init_scalar!(J_cpu, sim.mesh, params[:J])
    J = kernel_array(J_cpu)

    field = create_zeros(3 * n_total)

    ft = haskey(params, :ft) ? params[:ft] : t -> 1.0

    torque = SlonczewskiTorque(T(p[1]), T(p[2]), T(p[3]), T(beta), T(Lambda), T(xi), T(P),
                               J, field, ft, name)

    push!(sim.interactions, torque)

    return torque
end
