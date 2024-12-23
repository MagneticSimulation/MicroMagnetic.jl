export add_zeeman, update_zeeman, add_anis, update_anis, add_cubic_anis, add_exch, add_dmi,
       add_demag, add_dmi_int, add_exch_int, add_thermal_noise

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

Add a time varying zeeman to system.

The input `ft` is a function of time `t` and its return value should be a tuple with length 3.

Example:

```julia
  function time_fun(t)
    w = 2*pi*2.0e9
    return (sin(w*t), cos(w*t), 0)
  end

  function spatial_H(i, j, k, dx, dy, dz)
    H = 1e3
    if i<=2
        return (H, H, 0)
    end
    return (0, 0, 0)
  end

  add_zeeman(sim, spatial_H, time_fun)
```
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
    Hout = zeros(T, 3)

    zee = TimeZeeman(T(0), T(0), T(0), ft, field_kb, field, energy, name)
    push!(sim.interactions, zee)

    if sim.save_data
        id = length(sim.interactions)
        if isa(H0, Tuple) && length(H0) == 3  # FIXME: the output should depends on time!!!
            field_item = SaverItem((string(name, "_Hx"), string(name, "_Hy"),
                                    string(name, "_Hz")), ("<A/m>", "<A/m>", "<A/m>"),
                                   o::AbstractSim -> (H0[1]*zee.time_fx, H0[2]*zee.time_fy, H0[3]*zee.time_fz))
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
    add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
    E_\\mathrm{anis} =  K_{u} (1 - \\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0, 0, 1),
                  name="anis")
    n_total = sim.n_total
    T = Float[]
    Kus = zeros(T, n_total)
    init_scalar!(Kus, sim.mesh, Ku)

    Kus_kb = KernelAbstractions.zeros(default_backend[], T, n_total)
    field = KernelAbstractions.zeros(default_backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(default_backend[], T, n_total)

    lt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)
    naxis = (T(axis[1] / lt), T(axis[2] / lt), T(axis[3] / lt))

    copyto!(Kus_kb, Kus)
    anis = Anisotropy(Kus_kb, naxis, field, energy, name)
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
    add_cubic_anis(sim::AbstractSim, Kc::Float64; axis1=(1,0,0), axis2=(0,1,0), name="cubic")

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
    add_thermal_noise(sim::AbstractSim, Temp::NumberOrArrayOrFunction; name="thermal", scaling=t -> 1.0, k_B=k_B)

Adds thermal noise fields to the simulation. For micromagnetic model, the thermal noise is defined as

```math
\mathbf{b}^u = \eta \sqrt \frac{2 \alpha k_B T}{\mu_0 M_s \gamma \Delta V dt}
```
and $\eta$ is a random number follows the normal distribution where $\gamma=2.211\times 10^5$ m/(A·s) is the gyromagnetic ratio.

For the atomistic model, the thermal noise is defined as

```math
\mathbf{b}^u = \eta \sqrt \frac{2 \alpha k_B T}{\gamma \mu_s dt}.
```
where $\eta$ is a random number follows the normal distribution and $\gamma=1.76\times 10^{11}$ rad/(T·s). 

### Arguments
- `sim::AbstractSim`: The simulation object.
- `Temp::NumberOrArrayOrFunction`: Temperature value (can be a constant, array, or function).
- `name::String`: Name for the noise field (default: `"thermal"`).
- `scaling::Function`: A function to scale the noise over time (default: `t -> 1.0`).
- `k_B::Float64`: Boltzmann constant (default: `k_B`).

### Example
```julia
# Add thermal noise with a constant temperature of 100 K and a scaling function
add_thermal_noise(sim, 100.0, scaling = t -> exp(-t/1e-9))
```
"""
function add_thermal_noise(sim::AbstractSim, Temp::NumberOrArrayOrFunction; name="thermal",
                           scaling=t -> 1.0, k_B=k_B)
    N = sim.n_total
    T = Float[]
    field = KernelAbstractions.zeros(default_backend[], T, 3 * N)
    energy = KernelAbstractions.zeros(default_backend[], T, N)

    Spatial_T = KernelAbstractions.zeros(default_backend[], T, N)
    eta = KernelAbstractions.zeros(default_backend[], T, 3 * N)

    init_scalar!(Spatial_T, sim.mesh, Temp)
    thermal = StochasticField(Spatial_T, eta, field, energy, -1, name, k_B, scaling,
                              average(Spatial_T), T(1))

    push!(sim.interactions, thermal)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "<J>",
                        o::AbstractSim -> sum(o.interactions[id].energy)))

        push!(sim.saver.items,
              SaverItem(string("T_", name), "<K>",
                        o::AbstractSim -> o.interactions[id].average_temperature *
                                          o.interactions[id].scaling_factor))
    end
    return thermal
end

@doc raw"""
    add_exch_int(sim::AbstractSim, J::Float64; k1=1, k2=-1, name="rkky")

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
function add_exch_int(sim::MicroSim, J::Float64; k1=1, k2=-1, name="exch_int")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    if k2 == -1
        k2 = sim.mesh.nz
    end
    exch = InterlayerExchange(J, Int32(k1), Int32(k2), field, energy, name)

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
