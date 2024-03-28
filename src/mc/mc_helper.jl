abstract type MonteCarloEnergy end

mutable struct ExchangeMC{TF<:AbstractFloat} <: MonteCarloEnergy
    J::CuArray{TF, 1} #nearest exchange constant
    J1::CuArray{TF, 1} #next-nearest exchange constant
end

mutable struct NearestExchangeMC{TF<:AbstractFloat} <: MonteCarloEnergy
    J::CuArray{TF, 1} #nearest exchange constant
end

mutable struct Nearest_DMI_MC{TF<:AbstractFloat} <: MonteCarloEnergy
    D::CuArray{TF, 2} #nearest DMI vector
end

mutable struct ZeemanMC{TF<:AbstractFloat}
    Hx::TF
    Hy::TF
    Hz::TF
end

mutable struct DMI_MC{TF<:AbstractFloat} <: MonteCarloEnergy
    D::CuArray{TF, 2} #nearest DMI vector
    D1::CuArray{TF, 2} #next-nearest DMI vector
end

mutable struct UniformAnisotropyMC{TF<:AbstractFloat} <: MonteCarloEnergy
    Ku::TF
    ux::TF
    uy::TF
    uz::TF
    Kc::TF
end

mutable struct KagomeAnisotropyMC{TF<:AbstractFloat} <: MonteCarloEnergy
    Ku::TF
end


mutable struct MonteCarlo{TF<:AbstractFloat} <:AbstractSimGPU
  mesh::Mesh
  shape::CuArray{Bool, 1}
  exch::MonteCarloEnergy
  dmi::MonteCarloEnergy
  zee::ZeemanMC
  anis::MonteCarloEnergy
  saver::DataSaver
  spin::CuArray{TF, 1}
  nextspin::CuArray{TF, 1}
  rnd::CuArray{TF, 1}
  energy::CuArray{TF, 1}
  delta_E::CuArray{TF, 1}
  total_energy::TF
  n_nodes::Int64
  steps::Int64
  name::String
  T::Float64
  mc_2d::Bool
  MonteCarlo{T}() where {T<:AbstractFloat} = new()
end

function uniform_random_sphere(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    rand!(rnd)
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr uniform_random_sphere_kernel!(spin, rnd, N)

    return  nothing
end

function uniform_random_circle_xy(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    rand!(rnd)
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr uniform_random_circle_xy_kernel!(spin, rnd, N)

    return  nothing
end

function compute_site_energy_exch_dmi(sim::MonteCarlo, energy::ExchangeMC, bias::Int64, mesh::CubicMeshGPU)
    exch = sim.exch
    dmi = sim.dmi
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr add_dE_exch_dmi_cubic_mesh_kernel!(sim.spin, sim.nextspin,
                                        sim.shape,
                                        sim.energy, mesh.ngbs, mesh.nngbs,
                                        mesh.n_ngbs,
                                        exch.J, exch.J1, dmi.D, dmi.D1,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_exch_dmi(sim::MonteCarlo, energy::NearestExchangeMC, bias::Int64, mesh::TriangularMeshGPU)
    exch = sim.exch
    dmi = sim.dmi
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr add_dE_exch_dmi_triangular_mesh_kernel!(sim.spin, sim.nextspin,
                                        sim.shape,
                                        sim.energy, mesh.ngbs,
                                        mesh.n_ngbs,
                                        exch.J, dmi.D,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end


function compute_site_energy_zeeman_anis(sim::MonteCarlo, energy::UniformAnisotropyMC, bias::Int64, mesh::CubicMeshGPU)
    zee = sim.zee
    anis = sim.anis
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr dE_zeeman_anisotropy_cubic_mesh_kernel!(sim.spin,
                                        sim.nextspin, sim.shape,
                                        sim.energy, mesh.ngbs,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku, anis.Kc, anis.ux, anis.uy, anis.uz,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_zeeman_anis(sim::MonteCarlo, energy::UniformAnisotropyMC, bias::Int64, mesh::TriangularMeshGPU)
    zee = sim.zee
    anis = sim.anis
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr dE_zeeman_anisotropy_triangular_mesh_kernel!(sim.spin,
                                        sim.nextspin, sim.shape,
                                        sim.energy, mesh.ngbs,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku, anis.ux, anis.uy, anis.uz,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_zeeman_anis(sim::MonteCarlo, energy::KagomeAnisotropyMC, bias::Int64, mesh::TriangularMeshGPU)
    zee = sim.zee
    anis = sim.anis
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr dE_zeeman_kagome_anisotropy_triangular_mesh_kernel!(sim.spin,
                                        sim.nextspin, sim.shape,
                                        sim.energy, mesh.ngbs,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_single(sim::MonteCarlo, bias::Int64, mesh::AtomicMeshGPU)

    compute_site_energy_zeeman_anis(sim, sim.anis, bias, mesh) #include external fields as well
    compute_site_energy_exch_dmi(sim, sim.exch, bias, mesh) #include DMI as well

    return nothing
end

function compute_system_energy(sim::MonteCarlo, energy::UniformAnisotropyMC, mesh::AtomicMeshGPU)

    zee = sim.zee
    anis = sim.anis
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr total_E_zeeman_anisotropy_kernel!(sim.spin,
                                        sim.shape,
                                        sim.energy,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku, anis.Kc, anis.ux, anis.uy, anis.uz,
                                        mesh.n_nodes)

    return nothing
end


function compute_system_energy(sim::MonteCarlo, energy::ExchangeMC, mesh::CubicMeshGPU)

    exch = sim.exch
    dmi = sim.dmi
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr add_total_E_exch_dmi_cubic_mesh_kernel!(sim.spin,
                                        sim.shape,
                                        sim.energy, mesh.ngbs, mesh.nngbs,
                                        mesh.n_ngbs,
                                        exch.J, exch.J1, dmi.D, dmi.D1,
                                        mesh.n_nodes)

    return nothing
end


function compute_system_energy(sim::MonteCarlo, energy::NearestExchangeMC, mesh::TriangularMeshGPU)

    exch = sim.exch
    dmi = sim.dmi
    blk, thr = cudims(sim.n_nodes)
    @cuda blocks=blk threads=thr add_total_E_exch_dmi_triangular_mesh_kernel!(sim.spin,
                                        sim.shape,
                                        sim.energy, mesh.ngbs,
                                        mesh.n_ngbs,
                                        exch.J,  dmi.D,
                                        mesh.n_nodes)

    return nothing
end

function compute_system_energy(sim::MonteCarlo)

    compute_system_energy(sim, sim.anis, sim.mesh)
    compute_system_energy(sim, sim.exch, sim.mesh)

    sim.total_energy = sum(sim.energy)*k_B

    return sim.total_energy
end
