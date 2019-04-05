function Sim(mesh::MeshGPU; driver="LLG", name="dyn")
  nxyz = mesh.nx*mesh.ny*mesh.nz
  spin = cuzeros(FloatGPU,3*nxyz)
  prespin = cuzeros(FloatGPU,3*nxyz)
  field = cuzeros(FloatGPU,3*nxyz)
  energy = cuzeros(FloatGPU,nxyz)
  Ms = zeros(FloatGPU,nxyz)
  Ms_inv = cuzeros(FloatGPU,nxyz)
  driver = create_driver_gpu(driver, nxyz)

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> o.total_energy, average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
  interactions = []

  return MicroSimGPU(mesh, driver, saver, spin, prespin, field, energy, Ms, Ms_inv, nxyz, name, interactions)

end

function set_Ms(sim::MicroSim, fun_Ms::Function)
    mesh = sim.mesh
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        sim.Ms[id] = fun_Ms(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    return true
end

function set_Ms(sim::MicroSim, Ms::Number)
    for i =1:sim.nxyz
        sim.Ms[i] = Ms
    end
    return true
end

function init_m0(sim::AbstractSim, m0::Any; norm=true)
  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
    normalise(sim.prespin, sim.nxyz)
  end
  sim.spin[:] .= sim.prespin[:]
end
