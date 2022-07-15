mutable struct MicroSimFEM <: AbstractSim
    mesh::FEMesh
    driver::Driver
    saver::DataSaver
    spin::Array{Float64, 1}
    prespin::Array{Float64, 1}
    field::Array{Float64, 1}
    energy::Array{Float64, 1}
    Ms::Array{Float64, 1}
    pins::Array{Bool, 1}
    n_nodes::Int64
    n_cells::Int64
    name::String
    interactions::Array
    save_data::Bool
    MicroSimFEM() = new()
  end


  function Sim(mesh::FEMesh; driver="LLG", name="dyn", integrator="DormandPrince", save_data=true)

    sim = MicroSimFEM()

    sim.name = name
    sim.mesh = mesh
    n_nodes = mesh.number_nodes
    n_cells = mesh.number_cells
    sim.n_nodes = n_nodes
    sim.n_cells = n_cells
    sim.spin = zeros(Float64,3*n_nodes)
    sim.prespin = zeros(Float64,3*n_nodes)
    sim.field = zeros(Float64,3*n_nodes)
    sim.energy = zeros(Float64, n_nodes)

    sim.Ms = zeros(Float64, n_cells)
    sim.pins = zeros(Bool, n_nodes)
    sim.save_data = save_data

    if save_data
        headers = ["step", "E_total", ("m_x", "m_y", "m_z")]
        units = ["<>", "<J>",("<>", "<>", "<>")]
        results = [o::AbstractSim -> o.saver.nsteps,
                o::AbstractSim -> sum(o.energy),
                average_m]
        sim.saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
    end

    if driver!="none"
      if driver in ("LLG", "LLG_STT", "LLG_STT_CPP") && save_data
          saver = sim.saver
          insert!(saver.headers, 2, "time")
          insert!(saver.units, 2, "<s>")
          insert!(saver.results, 2, o::AbstractSim -> o.saver.t)
      end
      sim.driver = create_driver(driver, integrator, n_nodes)
    end

   sim.interactions = []
   #println(sim)
   return sim
end


function set_Ms(sim::MicroSimFEM, Ms::Number)
  sim.Ms .= Ms
  return true
end



function init_m0(sim::MicroSimFEM, m0::TupleOrArrayOrFunction; norm=true)

  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
      normalise(sim.prespin, sim.n_nodes)
  end

  if any(isnan, sim.prespin)
    error("NaN is given by the input m0!")
  end

  sim.spin[:] .= sim.prespin[:]
end
