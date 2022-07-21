mutable struct MicroSimFEM <: AbstractSim
    mesh::FEMesh
    driver::Driver
    saver::DataSaver
    spin::Array{Float64, 1}
    prespin::Array{Float64, 1}
    field::Array{Float64, 1}
    energy::Array{Float64, 1}
    Ms::Array{Float64, 1}
    L_mu::Array{Float64, 1}
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
    sim.L_mu = zeros(Float64, 3*n_nodes)
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
  compute_L_Ms!(sim.L_mu, sim.mesh, sim.Ms)
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


function add_demag(sim::MicroSimFEM; name="demag")
  demag = init_demag(sim)
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


"""
    add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::MicroSimFEM, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")
  
  Kus =  zeros(Float64, sim.n_cells)
  init_scalar!(Kus, sim.mesh, Ku)
  field = zeros(Float64, 3*sim.n_nodes)
  energy = zeros(Float64, sim.n_nodes)

  lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
  axes = zeros(Float64, 3*sim.n_cells)
  for i = 1:sim.n_cells
    axes[3*(i-1)+1] = axis[1]/lt
    axes[3*(i-1)+2] = axis[2]/lt
    axes[3*(i-1)+3] = axis[3]/lt
  end

  K_matrix = spzeros(3*sim.n_nodes, 3*sim.n_nodes)

  anis =  AnisotropyFEM(Kus, axes, field, energy, K_matrix, false, name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return anis
end


function save_inp(sim::MicroSimFEM, fname::String)
  mesh = sim.mesh

  f = open(fname, "w")

  data_number = 1 # FIX ME
  write(f, @sprintf("%d %d %d 0 0\n",
                   mesh.number_nodes,
                   mesh.number_cells,
                   data_number))# data number
  for n = 1:mesh.number_nodes
     write(f, @sprintf("%d %0.12g %0.12g %0.12g\n",n,
                      mesh.coordinates[1,n],
                      mesh.coordinates[2,n],
                      mesh.coordinates[3,n]))
  end

  for c = 1:mesh.number_cells
     write(f, @sprintf("%d %d tet %d %d %d %d\n",c,
                      1, #FIX ME: material id
                      mesh.cell_verts[1, c],
                      mesh.cell_verts[2, c],
                      mesh.cell_verts[3, c],
                      mesh.cell_verts[4, c]))
  end

  write(f, @sprintf("%d", data_number))
  for i=1:data_number
    write(f, " 3")
  end
  write(f, "\nm, none")

  m = reshape(sim.spin, 3, mesh.number_nodes)
  # FIX ME: to add more data
  for n = 1:mesh.number_nodes
    write(f, @sprintf("\n%d %0.12g %0.12g %0.12g",n, m[1,n], m[2,n], m[3,n]))
  end

 close(f)

end