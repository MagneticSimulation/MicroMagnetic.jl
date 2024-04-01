"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Add a static Zeeman energy to the simulation.
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    n_total = sim.n_total   
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    init_vector!(field, sim.mesh, H0)

    zeeman =  Zeeman(field, energy, name)
    push!(sim.interactions, zeeman)

    if sim.save_data
        id = length(sim.interactions)
        if isa(H0, Tuple)
            field_item = SaverItem(
                (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")),
                ("A/m", "A/m", "A/m"),
                o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
                )
            push!(sim.saver.items, field_item)
        end
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end

  @info "Standard Zeeman (CPU) has been added."

  return zeeman
end

"""
    update_zeeman(sim::AbstractSim, H0::Tuple; name="zeeman")

Set the Zeeman field to H0 where H0 is TupleOrArrayOrFunction according to its name. For example,

```julia
   add_zeeman(sim, (0,0,0), name="my_H")  #create a zeeman energy with field (0,0,0) A/m
   update_zeeman(sim, (0,0,1e5), name="my_H")  #change the field to (0,0,1e5) A/m
```

"""
function update_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    N_spins = sim.n_total
    field = zeros(Float64, 3*N_spins)
    init_vector!(field, sim.mesh, H0)

    for i in sim.interactions
        if i.name == name
            i.field[:] = field[:]
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
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")
    n_total = sim.n_total
    init_field = zeros(Float64, 3*n_total)
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    init_vector!(init_field, sim.mesh, H0)

    zeeman =  TimeZeeman(ft, init_field, field, energy, name)
    push!(sim.interactions, zeeman)

    if sim.save_data
        id = length(sim.interactions)
        if isa(H0, Tuple)
            field_item = SaverItem(
                (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")),
                ("<A/m>", "<A/m>", "<A/m>"),
                o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
                )
            push!(sim.saver.items, field_item)
        end
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return zeeman
end

"""
    add_exch_vector(sim::AbstractSim, A::TupleOrArrayOrFunction; name="exch")

Add a vector form exchange energy to the system. The exchange constant of 3 directions can be different.
For example:
```julia
add_exc_vector(sim, (2e-12,5e-12,0))
```
"""
function add_exch_vector(sim::AbstractSim, A::TupleOrArrayOrFunction; name="exch_vector")
    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    Spatial_A = zeros(Float64, 3*n_total)
    init_vector!(Spatial_A , sim.mesh, A)
    exch = Vector_Exchange(Spatial_A , field, energy, name)
    push!(sim.interactions, exch)
    
    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end

    return exch
end

"""
    add_exch(sim::AbstractSim, A::NumberOrArrayOrFunction; name="exch")

Add exchange energy to the system.
"""
function add_exch(sim::AbstractSim, A::NumberOrArrayOrFunction; name="exch")
    
    n_total = sim.n_total
    Spatial_A = zeros(Float64, sim.n_total)

  
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)

    init_scalar!(Spatial_A , sim.mesh, A)
    if isa(sim, MicroSim) 
        exch = Exchange(Spatial_A , field, energy, name)
    else
        exch = HeisenbergExchange(A, field, energy, name)
    end
    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return exch
end

""" 
    add_exch(sim::AbstractSim, geo::Geometry, A::Number; name="exch")

Add exchange interaction within the Geometry, or update corresponding A if other exch is added.
"""
function add_exch(sim::AbstractSim, geo::Geometry, A::Number; name="exch")
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.A, geo, A)
            return nothing
        end
    end
    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    Spatial_A = zeros(Float64, n_total)
    update_scalar_geometry(Spatial_A, geo, A)
    if isa(sim, MicroSim)
        exch = Exchange(Spatial_A , field, energy, name)
    else
        exch = HeisenbergExchange(A, field, energy, name)
    end
    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return exch
end

@doc raw"""
    add_exch_rkky(sim::AbstractSim, J::Float64; name="rkky")

Add an RKKY-type exchange to the system. The energy of RKKY-type exchange is defined as 

```math
E_\mathrm{rkky} =  - \int_\Gamma J_\mathrm{rkky} \mathbf{m}_{i} \cdot \mathbf{m}_{j} dA
```
where $\Gamma$ is the interface between two layers with magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$,
$J_\mathrm{rkky}$ is the coupling constant which is related to the spacer layer thickness. 

The effective field is given then as
```math
\mathbf{H}_i = \frac{1}{\mu_0 M_s}  \frac{J_\mathrm{rkky}}{\Delta} \mathbf{m}_{j} 
```
"""
function add_exch_rkky(sim::AbstractSim, J::Float64; name="rkky")
    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    exch = ExchangeRKKY(J, field, energy, name)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return exch
end

"""
    add_dmi(sim::AbstractSim, D::Tuple{Real, Real, Real}; name="dmi")

Add DMI to the system. Example:

```julia
   add_dmi(sim, (1e-3, 1e-3, 0))
```
"""
function add_dmi(sim::AbstractSim, D::Tuple{Real, Real, Real}; name="dmi")
    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    dmi =  BulkDMI(Float64(D[1]), Float64(D[2]), Float64(D[3]), field, energy, name)
    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return dmi
end

"""
    add_dmi(sim::AbstractSim, D::Real; name="dmi", type="bulk")

Add DMI to the system. `type` could be "bulk" or "interfacial"
Examples:

```julia
   add_dmi(sim, 1e-3, type="interfacial")
```
or
```julia
   add_dmi(sim, 1e-3, type="bulk")
```
"""
function add_dmi(sim::AbstractSim, D::Real; name="dmi", type="bulk")
    if type == "interfacial"
        return add_dmi_interfacial(sim, D, name=name)
    elseif type == "bulk"
        return add_dmi(sim, (D,D,D), name=name)
    else
        error("Supported DMI type:", "interfacial", "bulk")
    end

end

function add_dmi_interfacial(sim::AbstractSim, D::Real; name="dmi")
    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    dmi =  InterfacialDMI(Float64(D), field, energy, name)
    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return dmi
end

"""
    add_dmi(sim::AbstractSim;  name="dmi")

Add DMI to the system. Example:

```julia
   add_dmi(sim, (1e-3, 1e-3, 0))
```
"""
function add_dmi(sim::AbstractSim, name="dmi")
    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    z = 0.0
    dmi =  DMI(z,z,z,z,z,z,z,z,z,z,z,z,z,z, field, energy, name)
    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return dmi
end

"""
    add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)

Add Demag to the system. `Nx`, `Ny` and `Nz` can be used to describe the macro boundary conditions which means that
the given mesh is repeated `2Nx+1`, `2Ny+1 and `2Nz+1` times in `x`, `y` and `z` direction, respectively.
"""
function add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)
    demag = init_demag(sim, Nx, Ny, Nz)
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return demag
end

"""
    add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")
    n_total = sim.n_total
    Kus =  zeros(Float64, n_total)
    init_scalar!(Kus, sim.mesh, Ku)
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
    naxis = (axis[1]/lt, axis[2]/lt, axis[3]/lt)
    anis =  Anisotropy(Kus, naxis, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return anis
end

"""
    add_anis(sim::AbstractSim, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")

Add Anisotropy within the Geometry, or update corresponding Ku if other anis is added.
"""
function add_anis(sim::AbstractSim, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")
    lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
    naxis = (axis[1]/lt, axis[2]/lt, axis[3]/lt)
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.Ku, geo, Ku)
            return nothing
        end
    end
    n_total = sim.n_total
    Kus =  zeros(Float64, n_total)
    update_scalar_geometry(Kus, geo, Ku)
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    anis =  Anisotropy(Kus, naxis, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
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
function update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name = "anis")
    n_total = sim.n_total
    Kus =  zeros(Float64, n_total)
    init_scalar!(Kus, sim.mesh, Ku)
    for i in sim.interactions
        if i.name == name
            i.Ku[:] = Kus[:]
            return nothing
        end
    end
    return nothing
end

"""
    add_cubic_anis(sim::AbstractSim, Kc::Float64; axis1::Any=nothing, axis2::Any=nothing, name="cubic")

add cubic anisotropy with default axis (1,0,0) , (0,1,0), and (0,0,1)
use axis1=(1,1,0), axis2=(1,-1,0) to set a pair of normal axis, and the third axis will be calculated automatically
"""
function add_cubic_anis(sim::AbstractSim, Kc::Float64; axis1::Any=nothing, axis2::Any=nothing, name="cubic")
    axis1 = axis1 == nothing ? (1,0,0) : axis1
    axis2 = axis2 == nothing ? (0,1,0) : axis2
    norm1 = sqrt(axis1[1]^2+axis1[2]^2+axis1[3]^2)
    norm2 = sqrt(axis2[1]^2+axis2[2]^2+axis2[3]^2)
    naxis1,naxis2 = axis1./norm1,axis2./norm2
    if abs.(sum(naxis1.*naxis2)) > 1e-10
        @error("cubic axis not normal!")
        return nothing
    end
    naxis3= cross_product(axis1,axis2)
    axis = zeros(Float64,9)
    for i = 1:3
        axis[i] = naxis1[i]
        axis[i+3] = naxis2[i]
        axis[i+6] = naxis3[i]
    end

    n_total = sim.n_total
    field = zeros(Float64, 3*n_total)
    energy = zeros(Float64, n_total)
    anis =  CubicAnisotropy(axis, Kc, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSim->sum(o.interactions[id].energy)))
    end
    return anis
end