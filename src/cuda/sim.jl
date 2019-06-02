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


function init_m0(sim::MicroSimGPU, m0::Any; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end
  copyto!(sim.spin, spin)
  copyto!(sim.prespin, sim.spin)
  return true
end

function add_zeeman(sim::MicroSimGPU, H0::Any; name="zeeman")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  init_vector!(field, sim.mesh, H0)
  b = reshape(field, 3, nxyz)
  Hx = sum(b[1,:])/nxyz
  Hy = sum(b[2,:])/nxyz
  Hz = sum(b[3,:])/nxyz
  field_gpu = CuArray(field)
  zeeman =  ZeemanGPU(Float64(Hx), Float64(Hy), Float64(Hz), field, field_gpu, energy, Float(0.0), name)
  push!(sim.interactions, zeeman)

  push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
  push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
  id = length(sim.interactions)
  fun = o::AbstractSim ->  (o.interactions[id].Hx, o.interactions[id].Hy, o.interactions[id].Hz)
  push!(sim.saver.results, fun)

  push!(sim.saver.headers, string("E_",name))
  push!(sim.saver.units, "J")
  push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
  return zeeman
end


function add_exch(sim::MicroSimGPU, A::Float64; name="exch")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  exch = ExchangeGPU(Float(A), field, energy, Float(0.0), name)

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

function add_anis(sim::MicroSimGPU, Ku::Any; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Float = _cuda_using_double.x ? Float64 : Float32
  Kus = zeros(Float, sim.nxyz)
  init_scalar!(Kus, sim.mesh, Ku)
  Kus =  CuArray(Kus)
  field = zeros(Float, 3*nxyz)
  energy = zeros(Float, nxyz)
  anis =  AnisotropyGPU(Kus, axis, field, energy, Float(0.0), name)
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
