"""
    set_Ms_cylindrical(sim::MicroSimGPU, Ms::Number; axis=ez, r1=0, r2=0)
"""
function set_Ms_cylindrical(sim::MicroSimGPU, Ms::Number; axis=ez, r1=0, r2=0)
    geo = create_cylinder(sim.mesh, axis, r1=r1, r2=r2)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms_array = zeros(Float, sim.n_total)
    for i in 1:sim.n_total
        if geo.shape[i]
            Ms_array[i] = Ms
        end
    end
    copyto!(sim.Ms, Ms_array)
    return true
end


function set_pinning(sim::MicroSimGPU, ids::ArrayOrFunction)
    pins = zeros(Bool, sim.n_total)
    init_scalar!(pins, sim.mesh, ids)
    copyto!(sim.pins, pins)
    return true
end

function add_exch(sim::MicroSimGPU, A::NumberOrArrayOrFunction; name="exch")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    Spatial_A = CUDA.zeros(Float, n_total)
    init_scalar!(Spatial_A , sim.mesh, A)
    exch = ExchangeGPU(Spatial_A, field, energy, Float(0.0), name)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return exch
end

"""
    add_exch_anis(sim::MicroSimGPU, kea::NumberOrArrayOrFunction; name="exch_anis")

Add exchange anistropy to the system.
Ref: 10.1103/PhysRevResearch.2.043386
"""
function add_exch_anis(sim::MicroSimGPU, kea::NumberOrArrayOrFunction; name="exch_anis")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    Spatial_kea = CUDA.zeros(Float, n_total)
    init_scalar!(Spatial_kea , sim.mesh, kea)
    exch = ExchangeAnistropyGPU(Spatial_kea, field, energy, Float(0.0), name)
    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return exch
end

function add_exch(sim::MicroSimGPU, geo::Geometry, A::Number; name="exch")
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.A, geo, A)
            return nothing
        end
    end
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    Spatial_A = CUDA.zeros(Float, n_total)
    update_scalar_geometry(Spatial_A , geo, A)
    exch = ExchangeGPU(Spatial_A, field, energy, Float(0.0), name)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return exch
end


function add_thermal_noise(sim::MicroSimGPU, T::NumberOrArrayOrFunction; name="thermal", k_B=k_B)
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    Spatial_T = CUDA.zeros(Float, n_total)
    eta = CUDA.zeros(Float, 3*n_total)
    init_scalar!(Spatial_T , sim.mesh, T)
    thermal = StochasticFieldGPU(Spatial_T, eta, field, energy, Float(0.0), -1, name, k_B)

    push!(sim.interactions, thermal)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return thermal
end

function add_exch_vector(sim::MicroSimGPU, A::TupleOrArrayOrFunction; name="exch_vector")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    Spatial_A = CUDA.zeros(Float, 3*n_total)
    init_vector!(Spatial_A , sim.mesh, A)
    exch = Vector_ExchangeGPU(Spatial_A , field, energy, Float(0.0), name)
    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return exch
end

"""
    add_exch_rkky(sim::MicroSimGPU, J::Float64; name="rkky")

Add RKKY exchange energy to the system.
"""
function add_exch_rkky(sim::MicroSimGPU, J::Float64; name="rkky")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    exch = ExchangeRKKYGPU(Float(J), field, energy, Float(0.0), name)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return exch
end


function add_dmi(sim::MicroSimGPU, D::Tuple{Real, Real, Real}; name="dmi")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    dmi = BulkDMIGPU(Float(D[1]), Float(D[2]), Float(D[3]), field, energy, Float(0.0), name)

    push!(sim.interactions, dmi)
    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return dmi
end

@doc raw"""
    add_dmi_interlayer(sim::MicroSimGPU, D::Tuple{Real, Real, Real}; name="dmi_int")

Add an interlayer DMI to the system. The energy of interlayer DMI is defined as 

```math
E_\mathrm{dmi-int} =  \int_\Gamma \mathbf{D} \cdot \left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right) dA
```
where $\Gamma$ is the interface between two layers with magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$. 
$\mathbf{D}$ is the effective DMI vector. 

The effective field is given
```math
\mathbf{H}_i = \frac{1}{\mu_0 M_s \Delta}  \mathbf{D} \times \mathbf{m}_{j} 
```

"""
function add_dmi_interlayer(sim::MicroSimGPU, D::Tuple{Real, Real, Real}; name="dmi_int")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    dmi = InterlayerDMIGPU(Float(D[1]), Float(D[2]), Float(D[3]), field, energy, Float(0.0), name)

    push!(sim.interactions, dmi)
    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return dmi
end

function add_dmi(sim::MicroSimGPU, Dfun::Function; name="dmi")
    n_total = sim.n_total
    T = _cuda_using_double.x ? Float64 : Float32
    Ds = zeros(T, sim.n_total)
    init_scalar!(Ds, sim.mesh, Dfun)
    Ds_gpu = CuArray(Ds)
    field = zeros(T, 3*n_total)
    energy = zeros(T, n_total)
    dmi = SpatialBulkDMIGPU(Ds_gpu, field, energy, T(0.0), name)

    push!(sim.interactions, dmi)
    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return dmi
end

function add_dmi(sim::MicroSimGPU, D::Real; name="dmi", type="bulk")
    if type == "interfacial"
        return add_dmi_interfacial(sim, D, name=name)
    end
   return add_dmi(sim, (D,D,D), name=name)
end


function add_dmi_interfacial(sim::MicroSimGPU, D::Real; name="dmi")
    n_total = sim.n_total
    T = _cuda_using_double.x ? Float64 : Float32
    field = zeros(T, 3*n_total)
    energy = zeros(T, n_total)
    dmi =  InterfacialDMIGPU(T(D), field, energy, T(0.0), name)
    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return dmi
end

function add_dmi_interfacial(sim::MicroSimGPU, Dfun::Function; name="dmi")
    n_total = sim.n_total
    T = _cuda_using_double.x ? Float64 : Float32
    Ds = zeros(T, sim.n_total)
    init_scalar!(Ds, sim.mesh, Dfun)
    Ds_gpu = CuArray(Ds)
    field = zeros(T, 3*n_total)
    energy = zeros(T, n_total)
    dmi =  SpatialInterfacialDMIGPU(Ds_gpu, field, energy, T(0.0), name)
    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
    return dmi
end


function add_anis(sim::AbstractSimGPU, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    Kus = zeros(Float, sim.n_total)
    init_scalar!(Kus, sim.mesh, Ku)
    Kus =  CuArray(Kus)
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
    naxis = (axis[1]/lt,axis[2]/lt,axis[3]/lt)
    anis = AnisotropyGPU(Kus, naxis, field, energy, Float(0.0), name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
  @info "AnisotropyGPU has been added to the simulation."
  return anis
end

function add_anis(sim::AbstractSimGPU, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")

    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.Ku, geo, Ku)
            return nothing
        end
    end
    n_total = sim.n_total
    Float = _cuda_using_double.x ? Float64 : Float32
    Kus = zeros(Float, sim.n_total)
    update_scalar_geometry(Kus, geo, Ku)
    Kus =  CuArray(Kus)
    field = zeros(Float, 3*n_total)
    energy = zeros(Float, n_total)
    lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
    naxis = (axis[1]/lt,axis[2]/lt,axis[3]/lt)
    anis = AnisotropyGPU(Kus, naxis, field, energy, Float(0.0), name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
  return anis
end

#="""
    add_cubic_anis(sim::AbstractSimGPU, Kc::Float64, axis1::Tuple, axis2::Tuple; name="cubic")

Add cubic anistropy with two perpendicular axis axis1, axis2. (Only for SimGPU)

axis1 is a Tuple with size (3). 

Make sure axis1.axis2 == 0,and the third axis will be calculated automatically.

For example:

```julia
add_cubic_anis(sim, 1e3, (1, 1, 0), (1, -1, 0))
```
will add a cubic anistropy with axis (1, 1, 0), (1, -1, 0) and (0, 0, 1)
"""=#
function add_cubic_anis(sim::AbstractSimGPU, Kc::Float64; axis1::Any=nothing, axis2::Any=nothing, name="cubic")
  axis1 = axis1 == nothing ? (1,0,0) : axis1
  axis2 = axis2 == nothing ? (0,1,0) : axis2

  norm1 = sqrt(axis1[1]^2+axis1[2]^2+axis1[3]^2)
  norm2 = sqrt(axis2[1]^2+axis2[2]^2+axis2[3]^2)
  naxis1,naxis2 = axis1./norm1,axis2./norm2
  if abs.(sum(naxis1.*naxis2)) > 1e-10
    println("cubic axis not normal!")
    return nothing
  end
  naxis3= cross_product(axis1,axis2)

  n_total = sim.n_total
  Float = _cuda_using_double.x ? Float64 : Float32
  axis = zeros(Float,9)
  for i = 1:3
    axis[i] = naxis1[i]
    axis[i+3] = naxis2[i]
    axis[i+6] = naxis3[i]
  end
  field = zeros(Float, 3*n_total)
  energy = zeros(Float, n_total)
  anis = CubicAnisotropyGPU(axis, Float(Kc), field, energy, Float(0.0), name)
  push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items, SaverItem(string("E_", name), "J", o::AbstractSimGPU->o.interactions[id].total_energy))
    end
  return anis
end


function update_anis(sim::AbstractSimGPU, Ku::NumberOrArrayOrFunction; name = "anis")
  n_total = sim.n_total
  Kus =  zeros(Float64, n_total)
  init_scalar!(Kus, sim.mesh, Ku)
  Kus =  CuArray(Kus)
  for i in sim.interactions
    if i.name == name
      i.Ku[:] = Kus[:]
      return nothing
    end
  end
  return nothing
end

function add_demag(sim::MicroSimGPU; name="demag", Nx=0, Ny=0, Nz=0)
    demag = init_demag_gpu(sim, Nx, Ny, Nz)
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::MicroSimGPU->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return demag
end

function add_demag_gpu(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)
    demag = init_demag_gpu_II(sim, Nx, Ny, Nz)
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::MicroSim->sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return demag
end



