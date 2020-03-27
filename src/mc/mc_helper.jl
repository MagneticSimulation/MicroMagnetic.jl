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


mutable struct MonteCarloNew{TF<:AbstractFloat} <:AbstractSimGPU
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
  nxyz::Int64
  steps::Int64
  name::String
  T::Float64
  mc_2d::Bool
  MonteCarloNew{T}() where {T<:AbstractFloat} = new()
end

function uniform_random_sphere(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    rand!(rnd)
    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr uniform_random_sphere_kernel!(spin, rnd, N)

    return  nothing
end

function uniform_random_circle_xy(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    rand!(rnd)
    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr uniform_random_circle_xy_kernel!(spin, rnd, N)

    return  nothing
end

function compute_site_energy_cubic_mesh(sim::MonteCarloNew, energy::ExchangeMC, bias::Int64)
    mesh = sim.mesh
    exch = sim.exch
    dmi = sim.dmi
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr add_dE_exch_dmi_cubic_mesh_kernel!(sim.spin, sim.nextspin,
                                        sim.shape,
                                        sim.energy, mesh.ngbs, mesh.nngbs,
                                        mesh.n_ngbs,
                                        exch.J, exch.J1, dmi.D, dmi.D1,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_triangular_mesh(sim::MonteCarloNew, energy::NearestExchangeMC, bias::Int64)
    mesh = sim.mesh
    exch = sim.exch
    dmi = sim.dmi
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr add_dE_exch_dmi_triangular_mesh_kernel!(sim.spin, sim.nextspin,
                                        sim.shape,
                                        sim.energy, mesh.ngbs,
                                        mesh.n_ngbs,
                                        exch.J, dmi.D,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end


function compute_site_energy_cubic_mesh(sim::MonteCarloNew, energy::UniformAnisotropyMC, bias::Int64)

    mesh = sim.mesh
    zee = sim.zee
    anis = sim.anis
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr dE_zeeman_anisotropy_cubic_mesh_kernel!(sim.spin,
                                        sim.nextspin, sim.shape,
                                        sim.energy, mesh.ngbs,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku, anis.Kc, anis.ux, anis.uy, anis.uz,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_triangular_mesh(sim::MonteCarloNew, energy::UniformAnisotropyMC, bias::Int64)
    mesh = sim.mesh
    zee = sim.zee
    anis = sim.anis
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr dE_zeeman_anisotropy_triangular_mesh_kernel!(sim.spin,
                                        sim.nextspin, sim.shape,
                                        sim.energy, mesh.ngbs,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku, anis.ux, anis.uy, anis.uz,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_triangular_mesh(sim::MonteCarloNew, energy::KagomeAnisotropyMC, bias::Int64)
    mesh = sim.mesh
    zee = sim.zee
    anis = sim.anis
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr dE_zeeman_kagome_anisotropy_triangular_mesh_kernel!(sim.spin,
                                        sim.nextspin, sim.shape,
                                        sim.energy, mesh.ngbs,
                                        zee.Hx, zee.Hy, zee.Hz,
                                        anis.Ku,
                                        mesh.nx, mesh.ny, mesh.nz, bias)
end

function compute_site_energy_single(sim::MonteCarloNew, bias::Int64, cubic::Bool)

    if cubic
        compute_site_energy_cubic_mesh(sim, sim.anis, bias) #include external fields as well
        compute_site_energy_cubic_mesh(sim, sim.exch, bias) #include DMI as well
    else
        compute_site_energy_triangular_mesh(sim, sim.anis, bias)  #include external fields as well
        compute_site_energy_triangular_mesh(sim, sim.exch, bias)  #include DMI as well
    end

    return nothing
end

function compute_system_energy(sim::MonteCarloNew)

    sim.nextspin .= sim.spin
    sim.spin .= 0

    compute_site_energy_single(sim, 0)
    compute_site_energy_single(sim, 1)
    compute_site_energy_single(sim, 2)

    sim.spin .= sim.nextspin

    return sum(sim.energy)*k_B
end
