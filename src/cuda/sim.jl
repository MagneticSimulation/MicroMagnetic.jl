function Sim(mesh::MeshGPU; driver="LLG", name="dyn")
  nxyz = mesh.nx*mesh.ny*mesh.nz
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = cuzeros(Float, 3*nxyz)
  prespin = cuzeros(Float,3*nxyz)
  field = cuzeros(Float,3*nxyz)
  energy = cuzeros(Float,nxyz)
  Ms = cuzeros(Float, nxyz)
  driver = create_driver_gpu(driver, nxyz)

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> o.total_energy, average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
  interactions = []

  blocks, threads = CuArrays.cudims(nxyz)
  return MicroSimGPU(mesh, driver, saver, spin, prespin, field, energy, Ms, Float(0.0),
                     nxyz, blocks, threads, name, interactions)

end

function set_Ms(sim::MicroSimGPU, init::Any)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms = zeros(Float, sim.nxyz)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.Ms, Ms)
    return true
end


function init_m0(sim::MicroSimGPU, m0::TupleOrArrayOrFunction; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end
  Ms = Array(sim.Ms)
  for i = 1:sim.nxyz
      if abs(Ms[i]) < eps(Float)
          spin[3*i-2] = 0
          spin[3*i-1] = 0
          spin[3*i] = 0
      end
  end
  copyto!(sim.spin, spin)
  copyto!(sim.prespin, sim.spin)
  return true
end

function add_zeeman(sim::MicroSimGPU, H0::TupleOrArrayOrFunction; name="zeeman")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  init_vector!(field, sim.mesh, H0)

  field_gpu = CuArray(field)
  zeeman =  ZeemanGPU(field, energy, field_gpu, T(0.0), name)
  push!(sim.interactions, zeeman)

  if isa(H0, Tuple)
      push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
      push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
      id = length(sim.interactions)
      fun = o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
      push!(sim.saver.results, fun)
  end

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return zeeman
end

function update_zeeman(z::ZeemanGPU, H0::Tuple)
  b = reshape(z.field, 3, div(length(z.field),3))
  b[1, :] .= H0[1]
  b[2, :] .= H0[2]
  b[3, :] .= H0[3]
  copyto!(z.cufield, z.field)
  return nothing
end

function add_zeeman(sim::MicroSimGPU, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  init_vector!(field, sim.mesh, H0)
  init_field = CuArray(field)

  zeeman =  TimeZeemanGPU(ft, init_field, field, energy, T(0), name)
  push!(sim.interactions, zeeman)

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
  return zeeman
end

function add_exch(sim::MicroSimGPU, A::NumberOrArrayOrFunction; name="exch")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  Spatial_A = cuzeros(Float, nxyz)
  init_scalar!(Spatial_A , sim.mesh, A)
  exch = ExchangeGPU(Spatial_A, field, energy, Float(0.0), name)

  push!(sim.interactions, exch)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return exch
end

function add_exch_rkky(sim::MicroSimGPU, sigma::Float64, Delta::Float64; name="rkky")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  exch = ExchangeRKKYGPU(Float(sigma), Float(Delta), field, energy, Float(0.0), name)

  push!(sim.interactions, exch)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return exch
end


function add_dmi(sim::MicroSimGPU, D::Tuple{Real, Real, Real}; name="dmi")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  dmi = BulkDMIGPU(Float(D[1]), Float(D[2]), Float(D[3]), field, energy, Float(0.0), name)

  push!(sim.interactions, dmi)
  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
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
  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
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

    push!(sim.saver.headers, string("E_",name))
    push!(sim.saver.units, "J")
    id = length(sim.interactions)
    push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    return dmi
end

function add_anis(sim::MicroSimGPU, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")
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

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return anis
end

function add_demag(sim::MicroSimGPU; name="demag")
  demag = init_demag_gpu(sim)
  demag.name = name
  push!(sim.interactions, demag)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::MicroSimGPU->o.interactions[id].total_energy)
  return demag
end

function add_demag_gpu(sim::MicroSim; name="demag")
  demag = init_demag_gpu_II(sim)
  demag.name = name
  push!(sim.interactions, demag)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  id = length(sim.interactions)
  push!(sim.saver.results, o::MicroSim->sum(o.interactions[id].energy))
  return demag
end
