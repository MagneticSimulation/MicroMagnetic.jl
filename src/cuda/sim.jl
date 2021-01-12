"""
    Sim(mesh::MeshGPU; driver="LLG", name="dyn", integrator="DormandPrince")

Create a simulation instance for GPU mesh.
"""
function Sim(mesh::MeshGPU; driver="LLG", name="dyn", integrator="DormandPrince", save_data=true)
    Float = _cuda_using_double.x ? Float64 : Float32
    if isa(mesh, FDMeshGPU)
        sim = MicroSimGPU{Float}()
    else
        sim = AtomicSimGPU{Float}()
    end
    sim.mesh = mesh
    nxyz = mesh.nx*mesh.ny*mesh.nz

    sim.nxyz = nxyz
    sim.spin = CUDA.zeros(Float, 3*nxyz)
    sim.prespin = CUDA.zeros(Float,3*nxyz)
    sim.field = CUDA.zeros(Float,3*nxyz)
    sim.energy = CUDA.zeros(Float,nxyz)
    if isa(mesh, FDMeshGPU)
        sim.Ms = CUDA.zeros(Float, nxyz)
    else
        sim.mu_s = CUDA.zeros(Float, nxyz)
    end
    sim.pins = CUDA.zeros(Bool, nxyz)
    sim.total_energy = 0.0
    sim.interactions = []
    sim.save_data = save_data

    if save_data
        headers = ["step", "E_total", ("m_x", "m_y", "m_z")]
        units = ["<>", "<J>",("<>", "<>", "<>")]
        results = [o::AbstractSim -> o.saver.nsteps,
                   o::AbstractSim -> o.total_energy, average_m]
        sim.saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
    end

    if driver!="none"
        sim.driver = create_driver_gpu(driver, integrator, nxyz)

        if driver in ("LLG", "LLG_STT", "LLG_STT_CPP") && save_data
            saver = sim.saver
            insert!(saver.headers, 2, "time")
            insert!(saver.units, 2, "<s>")
            insert!(saver.results, 2, o::AbstractSim -> o.saver.t)
        end
    end

    blocks, threads = cudims(nxyz)
    sim.blocks = blocks
    sim.threads = threads
    sim.name = name

    return sim

end

"""
    set_Ms(sim::MicroSimGPU, Ms::NumberOrArrayOrFunction)
"""
function set_Ms(sim::MicroSimGPU, Ms::NumberOrArrayOrFunction)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms_a = zeros(Float, sim.nxyz)
    init_scalar!(Ms_a, sim.mesh, Ms)
    copyto!(sim.Ms, Ms_a)
    return true
end


"""
    set_Ms_cylindrical(sim::MicroSimGPU, Ms::Number; axis=ez, r1=0, r2=0)
"""
function set_Ms_cylindrical(sim::MicroSimGPU, Ms::Number; axis=ez, r1=0, r2=0)
    geo = create_cylinder(sim.mesh, axis, r1=r1, r2=r2)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms_array = zeros(Float, sim.nxyz)
    for i in 1:sim.nxyz
        if geo.shape[i]
            Ms_array[i] = Ms
        end
    end
    copyto!(sim.Ms, Ms_array)
    return true
end


function set_pinning(sim::MicroSimGPU, ids::ArrayOrFunction)
    pins = zeros(Bool, sim.nxyz)
    init_scalar!(pins, sim.mesh, ids)
    copyto!(sim.pins, pins)
    return true
end

"""
    init_m0(sim::AbstractSimGPU, m0::TupleOrArrayOrFunction; norm=true)

Set the initial magnetization of the system. If `norm=false` the magnetization array will be not normalised.
"""
function init_m0(sim::AbstractSimGPU, m0::TupleOrArrayOrFunction; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end

  Ms = isa(sim.mesh, FDMeshGPU) ? Array(sim.Ms) : Array(sim.mu_s)
  for i = 1:sim.nxyz
      if Ms[i] == 0
          spin[3*i-2] = 0
          spin[3*i-1] = 0
          spin[3*i] = 0
      end
  end
  copyto!(sim.spin, spin)
  copyto!(sim.prespin, sim.spin)
  return true
end

"""
    add_zeeman(sim::AbstractSimGPU, H0::TupleOrArrayOrFunction; name="zeeman")

Add a static Zeeman energy to the simulation.
"""
function add_zeeman(sim::AbstractSimGPU, H0::TupleOrArrayOrFunction; name="zeeman")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  init_vector!(field, sim.mesh, H0)

  field_gpu = CuArray(field)
  zeeman =  ZeemanGPU(field, energy, field_gpu, T(0.0), name)
  push!(sim.interactions, zeeman)

  if sim.save_data
      if isa(H0, Tuple)
          push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
          push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
          id = length(sim.interactions)
          fun = o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
          push!(sim.saver.results, fun)
      end

      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return zeeman
end


function update_zeeman(sim::MicroSimGPU, H0::TupleOrArrayOrFunction; name = "zeeman")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  init_vector!(field, sim.mesh, H0)

  field_gpu = CuArray(field)

  for i in sim.interactions
    if i.name == name
      i.field[:] = field[:]
      i.cufield[:] = field_gpu[:]
      return nothing
    end
  end
  return nothing
end


function add_zeeman(sim::MicroSimGPU, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  local_filed = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  init_vector!(field, sim.mesh, H0)
  init_field = CuArray(field)
  cufield = CuArray(local_filed)

  zeeman =  TimeZeemanGPU(ft, init_field, cufield, field, energy, T(0), name)
  push!(sim.interactions, zeeman)

  if sim.save_data
      if isa(H0, Tuple)
          push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
          push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
          id = length(sim.interactions)
          fun = o::AbstractSim ->  (o.interactions[id].cufield[1], o.interactions[id].cufield[2], o.interactions[id].cufield[3])
          push!(sim.saver.results, fun)
      end

      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return zeeman
end

function add_exch(sim::MicroSimGPU, A::NumberOrArrayOrFunction; name="exch")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  Spatial_A = CUDA.zeros(Float, nxyz)
  init_scalar!(Spatial_A , sim.mesh, A)
  exch = ExchangeGPU(Spatial_A, field, energy, Float(0.0), name)

  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return exch
end

"""
    add_exch_anis(sim::MicroSimGPU, kea::NumberOrArrayOrFunction; name="exch_anis")

Add exchange anistropy to the system.
Ref: 10.1103/PhysRevResearch.2.043386
"""
function add_exch_anis(sim::MicroSimGPU, kea::NumberOrArrayOrFunction; name="exch_anis")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  Spatial_kea = CUDA.zeros(Float, nxyz)
  init_scalar!(Spatial_kea , sim.mesh, kea)
  exch = ExchangeAnistropyGPU(Spatial_kea, field, energy, Float(0.0), name)
  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
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
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  Spatial_A = CUDA.zeros(Float, nxyz)
  update_scalar_geometry(Spatial_A , geo, A)
  exch = ExchangeGPU(Spatial_A, field, energy, Float(0.0), name)

  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return exch
end


function add_thermal_noise(sim::MicroSimGPU, T::NumberOrArrayOrFunction; name="thermal")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  Spatial_T = CUDA.zeros(Float, nxyz)
  eta = CUDA.zeros(Float, 3*nxyz)
  init_scalar!(Spatial_T , sim.mesh, T)
  thermal = StochasticFieldGPU(Spatial_T, eta, field, energy, Float(0.0), -1, name)

  push!(sim.interactions, thermal)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return thermal
end

function add_exch_vector(sim::MicroSimGPU, A::TupleOrArrayOrFunction; name="exch_vector")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  Spatial_A = CUDA.zeros(Float, 3*nxyz)
  init_vector!(Spatial_A , sim.mesh, A)
  exch = Vector_ExchangeGPU(Spatial_A , field, energy, Float(0.0), name)
  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return exch
end

function add_exch_rkky(sim::MicroSimGPU, sigma::Float64, Delta::Float64; name="rkky")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  exch = ExchangeRKKYGPU(Float(sigma), Float(Delta), field, energy, Float(0.0), name)

  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return exch
end


function add_dmi(sim::MicroSimGPU, D::Tuple{Real, Real, Real}; name="dmi")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  dmi = BulkDMIGPU(Float(D[1]), Float(D[2]), Float(D[3]), field, energy, Float(0.0), name)

  push!(sim.interactions, dmi)
  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return dmi
end

function add_dmi(sim::MicroSimGPU, Dfun::Function; name="dmi")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  Ds = zeros(T, sim.nxyz)
  init_scalar!(Ds, sim.mesh, Dfun)
  Ds_gpu = CuArray(Ds)
  field = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  dmi = SpatialBulkDMIGPU(Ds_gpu, field, energy, T(0.0), name)

  push!(sim.interactions, dmi)
  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
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
    nxyz = sim.nxyz
    T = _cuda_using_double.x ? Float64 : Float32
    field = zeros(T, 3*nxyz)
    energy = zeros(T, nxyz)
    dmi =  InterfacialDMIGPU(T(D), field, energy, T(0.0), name)
    push!(sim.interactions, dmi)

    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
    return dmi
end

function add_dmi_interfacial(sim::MicroSimGPU, Dfun::Function; name="dmi")
    nxyz = sim.nxyz
    T = _cuda_using_double.x ? Float64 : Float32
    Ds = zeros(T, sim.nxyz)
    init_scalar!(Ds, sim.mesh, Dfun)
    Ds_gpu = CuArray(Ds)
    field = zeros(T, 3*nxyz)
    energy = zeros(T, nxyz)
    dmi =  SpatialInterfacialDMIGPU(Ds_gpu, field, energy, T(0.0), name)
    push!(sim.interactions, dmi)

    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
    return dmi
end


function add_anis(sim::AbstractSimGPU, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  Kus = zeros(Float, sim.nxyz)
  init_scalar!(Kus, sim.mesh, Ku)
  Kus =  CuArray(Kus)
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
  naxis = (axis[1]^2/lt,axis[2]^2/lt,axis[3]^2/lt)
  anis = AnisotropyGPU(Kus, naxis, field, energy, Float(0.0), name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return anis
end

function add_anis(sim::AbstractSimGPU, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")

  for interaction in sim.interactions
      if interaction.name == name
         update_scalar_geometry(interaction.Ku, geo, Ku)
         return nothing
      end
  end
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  Kus = zeros(Float, sim.nxyz)
  update_scalar_geometry(Kus, geo, Ku)
  Kus =  CuArray(Kus)
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
  naxis = (axis[1]^2/lt,axis[2]^2/lt,axis[3]^2/lt)
  anis = AnisotropyGPU(Kus, naxis, field, energy, Float(0.0), name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
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

  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  axis = zeros(Float,9)
  for i = 1:3
    axis[i] = naxis1[i]
    axis[i+3] = naxis2[i]
    axis[i+6] = naxis3[i]
  end
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  anis = CubicAnisotropyGPU(axis, Float(Kc), field, energy, Float(0.0), name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  end
  return anis
end


function update_anis(sim::AbstractSimGPU, Ku::NumberOrArrayOrFunction; name = "anis")
  nxyz = sim.nxyz
  Kus =  zeros(Float64, nxyz)
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
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::MicroSimGPU->o.interactions[id].total_energy)
  end
  return demag
end

function add_demag_gpu(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)
  demag = init_demag_gpu_II(sim, Nx, Ny, Nz)
  demag.name = name
  push!(sim.interactions, demag)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::MicroSim->sum(o.interactions[id].energy))
  end
  return demag
end
