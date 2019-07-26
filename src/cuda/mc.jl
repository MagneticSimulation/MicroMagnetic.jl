using Random

mutable struct MonteCarlo{TF<:AbstractFloat} <:AbstractSimGPU
  mesh::Mesh
  mu_s::Float64
  J::Float64
  J1::Float64
  D::Float64
  D1::Float64
  bulk_dmi::Bool
  Ku::Float64
  T::Float64
  Hx::Float64
  Hy::Float64
  Hz::Float64
  saver::DataSaver
  spin::CuArray{TF, 1}
  nextspin::CuArray{TF, 1}
  rnd::CuArray{TF, 1}
  energy::CuArray{TF, 1}
  total_energy::TF
  nxyz::Int64
  steps::Int64
  name::String
end

"""
Create a MonteCarlo simulation.

J, D and Ku in units of Joule
H in units of Tesla.

"""
function MonteCarlo(mesh::Mesh; mu_s=1.0*mu_B, J=1*meV, J1=0, D=0, D1=0, Ku=0, H=(0,0,0), T=10, bulk_dmi=false, name="dyn")
  nxy = mesh.nx*mesh.ny*mesh.nz
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = CuArrays.zeros(Float, 3*nxy)
  nextspin = CuArrays.zeros(Float,3*nxy)
  rnd = CuArrays.zeros(Float,3*nxy)
  energy = CuArrays.zeros(Float,nxy)
  if !(isa(mesh, TriangularMesh) || isa(mesh, CubicMeshGPU))
      error("MonteCarlo only support TriangularMesh or CubicMeshGPU, but input mesh is ", typeof(mesh))
  end

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> o.total_energy, average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)

  return MonteCarlo(mesh, mu_s, J/k_B, J1/k_B, D/k_B, D1/k_B, bulk_dmi, Ku/k_B, Float64(T),
                    H[1]*mu_s/k_B, H[2]*mu_s/k_B, H[3]*mu_s/k_B,
                    saver, spin, nextspin, rnd, energy, Float(0.0), nxy, 0, name)
end


function init_vector!(v::Array{T, 1}, mesh::TriangularMesh, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny
  dx,dy = mesh.dx, mesh.dy
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
        id = index(i, j, mesh.nx, mesh.ny)
        b[:, id] .=  init(i,j,dx,dy)
      end
  end
  return nothing
end

function init_vector!(v::Array{T, 1}, mesh::CubicMeshGPU, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
      for k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        b[:, id] .=  init(i,j,k)
      end
    end
  end
end

function init_m0(sim::MonteCarlo, m0::Any; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end
  copyto!(sim.spin, spin)
  return true
end


function uniform_random_sphere(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    rand!(rnd)
    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr uniform_random_sphere_kernel!(spin, rnd, N)

    return  nothing
end

function run_single_step_triangular(sim::MonteCarlo, bias::Int64)
  F = _cuda_using_double.x ? Float64 : Float32
  blk, thr = CuArrays.cudims(sim.nxyz)
  @cuda blocks=blk threads=thr run_step_triangular_kernel!(sim.spin, sim.nextspin,
                               sim.rnd, sim.energy, sim.mesh.ngbs,
                               F(sim.J), F(sim.J1), F(sim.D), F(sim.D1), sim.bulk_dmi, F(sim.Ku),
                               F(sim.Hx), F(sim.Hy), F(sim.Hz), F(sim.T),
                               sim.mesh.nx, sim.mesh.ny, bias)
  return nothing
end

function run_single_step_cubic(sim::MonteCarlo, bias::Int64)
  F = _cuda_using_double.x ? Float64 : Float32
  blk, thr = CuArrays.cudims(sim.nxyz)
  mesh = sim.mesh
  @cuda blocks=blk threads=thr run_step_cubic_kernel!(sim.spin, sim.nextspin,
                               sim.rnd, sim.energy, mesh.ngbs, mesh.nngbs,
                               F(sim.J), F(sim.J1), F(sim.D), F(sim.D1), sim.bulk_dmi, F(sim.Ku),
                               F(sim.Hx), F(sim.Hy), F(sim.Hz), F(sim.T),
                               mesh.nx, mesh.ny, mesh.nz, bias)
  return nothing
end


function run_step_triangular(sim::MonteCarlo)

  uniform_random_sphere(sim.nextspin, sim.rnd, sim.nxyz)
  run_single_step_triangular(sim, 0)
  run_single_step_triangular(sim, 1)
  run_single_step_triangular(sim, 2)
  sim.steps += 1

  return  nothing
end

function run_step_cubic(sim::MonteCarlo)

  uniform_random_sphere(sim.nextspin, sim.rnd, sim.nxyz)
  run_single_step_cubic(sim, 0)
  run_single_step_cubic(sim, 1)
  run_single_step_cubic(sim, 2)
  sim.steps += 1

  return  nothing
end

function compute_triangular_energy(sim::MonteCarlo)
    F = _cuda_using_double.x ? Float64 : Float32
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr total_energy_triangular_kernel!(sim.spin, sim.energy,
                                 sim.mesh.ngbs,
                                 F(sim.J), F(sim.J1), F(sim.D), F(sim.D1),
                                 sim.bulk_dmi, F(sim.Ku),
                                 F(sim.Hx), F(sim.Hy), F(sim.Hz),
                                 sim.mesh.nxyz)
    return sum(sim.energy)*k_B
end

function compute_cubic_energy(sim::MonteCarlo)
    F = _cuda_using_double.x ? Float64 : Float32
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr total_energy_cubic_kernel!(sim.spin, sim.energy,
                                 sim.mesh.ngbs, sim.mesh.nngbs,
                                 F(sim.J), F(sim.J1), F(sim.D), F(sim.D1),
                                 sim.bulk_dmi, F(sim.Ku),
                                 F(sim.Hx), F(sim.Hy), F(sim.Hz),
                                 sim.mesh.nxyz)
    return sum(sim.energy)*k_B
end

function run_sim(sim::MonteCarlo; maxsteps=10000, save_m_every = 10, save_vtk_every=-1)
    cubic = isa(sim.mesh, CubicMeshGPU) ? true : false
    for i=1:maxsteps

        if save_m_every>0
            if sim.steps%save_m_every == 0
                energy = cubic ? compute_cubic_energy(sim) : compute_triangular_energy(sim)
                @info @sprintf("step=%5d  total_energy=%g", sim.steps, energy)
                #compute_system_energy(sim, sim.spin, 0.0)
                #write_data(sim)
            end
        end

        if save_vtk_every > 0
            if sim.steps%save_vtk_every == 0
              save_vtk(sim, @sprintf("%s_%d", sim.name, sim.steps))
            end
        end

        if cubic
           run_step_cubic(sim)
        else
           run_step_triangular(sim)
        end

    end
end
